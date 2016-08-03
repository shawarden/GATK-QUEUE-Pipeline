/**
 * FastQ file to Gnomic VCF file.
 *
 * This scala script runs through GATK's Queue process to run alignment process in a multi-core and multi-node friendly manner.
 *
 * It takes a file containing a list of id requests, matches them to the sample description file line and passes these lines to the alignment through to variant calling.
 *
 * Command line options are:
 *  -a   or --alignemtns_file [file_path]: file containing alignment requests.
 *  -d   or --sample_description_file [file_path]: file containing sample file paths and known information such as capture platform, gender, pedigree, etc. Defaults settings file entry when not specified.
 *  -r   or --reference_file [file_path]: genomic reference file to use. Defaults settings file entry when not specified.
 *  -nt  or --number_of_threads [integer]: number of cpu threads to use. This is multiplicative for memory use as each thread has it's own memory (ram) block. Defaults to 4 when not specified.
 *  -nct or --number_of_core_threads [integer]: number of cpu cores to use. In theory this should allow processes with a single core to use the same memory (ram) block, reducing the overall memory consumed by the process. For now, this option is tied to -nt as most processes that can have -nt cannot use -nct and vise versa and having to specify both pointless.
 *  -nsc or --number_of_scatters [integer]: number of sub-jobs to allow processes to scatter to. Dividing up large jobs into smaller segments and running them in parallel instead of sequentially. Please note that only some functions can be scattered and the sum of all scattered jobs is greater than the sum of the single job as the data must be initially parcelled out to each process and then re-combined at the end of the process runs. Please also note that just because you've split a job up, does not mean that they will ALL run simultaneously. It will work on however many jobs it has been allowed to at a given instant. For now, using -nt will result in some degree of scatter-gather processing for functions that cannot use -nt or -nsc. Defaults to 1 when not specified and uses -nt number to generate reduced scatter count.
 *  -bwa or --bwa_path [file_path]: command path for Burrows Wheeler Aligner. Defaults settings file entry when not specified.
 *  -qm  or --qualimap_path [file_path]: command path for QualiMap. Defaults settings file entry when not specified.
 *  -xp  or --allow_cross_platform: no arguments given. this will allow the merging of Libraries from different exome capture platforms which is usually not permitted. When this is specified, quality-mapping will be done at the library level instead of individual level to give the best result. This will result in multiple quality-mapping runs. Defaults to not allowing cross capture platform merging at library level when this is not specified.
 *  -mp  or --mitochondira_ploidy [integer]: specified the ploidy at which to process mitochondiral data. Anything Less than 2 means excluded. Defaults to 0 if not specified.
 *  -s   or --file_folder_structure [folder|single]: specifies whether to have all your ouput in a single individual's folder or subdivided into dna, library and run folders to separate out the input and output files. Essensially replaces the folder / with a _ if anything other than 'folder' is set. Defaults to 'folder' if not specified.
 *
 * It takes various settings from a global settings file located at [/home/clinical_genetics/Resources/QSettings]
 * This includes:
 *   Genomic regions (PAR1, True X, Y, etc)
 *   Reference file for the human genome.
 *   Sample description file.
 *   Common variants file used for finger-printing runs and samples.
 *   Capture platform paths.
 *   Burrows Wheeler Aligner command path.
 *   QualiMap command path.
 *   Several pipeline control command paths.
 *   Mills, HapMap and DBSNP file paths.
 *   Additional known targets file paths which will ideally build over time.
 *
 * Prior to running the primary alignment phase (BWA alignment on the fastq file, Picard Sort Sam on that output and then Picard Clean Sam on the sorted output) it will verify sample description file data.
 *  It will:
 *    Check the settings file for any unspecified configurations required.
 *    Build a list of capture platforms and their gender-specific boundries.
 *    Separate jobs list by individual, library and runs.
 *    Attempt to locate trimmed files for the given sample.
 *    Verify sample desccription data all samples for an individual. (Gender/Capture platform consistenty)
 *
 * After Primary alignemnt of each run it will verify that all aligned, sorted and cleaned samples list as for a given individual are actually from the same individual. If the genotype concordance between any runs for an individual are inconsistent the pipeline for this individual will break and can be reviewed/restarted later once the samples have been verified, swapped or discarded.
 *
 * Once all runs are aligned and verified to be the same individual, the runs will be merged into a single library with GATK Realigner Target Creator, GATK Indel Realigner and finally Picard Mark Duplicates).
 *
 * If the capture platform is unknown, or there are multiple capture platforms for the various libraried for this individual, qualimap will run per library.
 *
 * If the capture platform is unknown for any sample or if cross platform merging is not permitted, the pipeline for this individual will break at this point to be resumed laster.
 *
 * If cross platform merging is permitted and all platforms are known, qualimap will still run per-library and they will continue on.
 *
 * 
 
 */
 
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._

import java.io._
//import scala.collection.JavaConversions._

class FastQtoGVCF extends org.broadinstitute.gatk.queue.QScript {
	
	qscript =>
	
	// File containing sequence(s) to be processed.
	// File should have 1 request per line with tab or space separated fields.
	@Input(
		doc="File containing alignment requests to lookup and process.",
		shortName="a",
		fullName="alignments_file",
		required=false
	)
	var alignmentsFile: File = "align.txt"
	
	@Argument(
		doc="Sample Descriptiosn file containing file paths",
		shortName="d",
		fullName="sample_description_file",
		required=false
	)
	var sampleDescriptionFile: File = _
	
	// Reference to use for primary alignment and indel realignment.
	// BWA doesn't like when it ends in .fastq so strip that off for BWA runs.
	@Argument(
		doc="FastA reference file to align against",
		shortName="r",
		fullName="reference_file",
		required=false
	)
	var referenceFile: File = _
	
	// Base output path for files to be created at.
	@Argument(
		doc="File folder location to output final data to. Some intermediate steps may occur on the system's local storage in your home folder/queueTMP",
		shortName="b",
		fullName="base_path",
		required=true
	)
	var basePath: File = _
	
	@Input(
		doc="Settings file location",
		shortName="c",
		fullName="config_file",
		required=true
	)
	var settingsFile: File = "settings.txt"
	
	// How many threads can the tools use?
	// This will be passed to BWA and GATK tools.
	@Argument(
		doc="Integer number of threads to allocate to processes.",
		shortName="nt",
		fullName="number_of_threads",
		required=false
	)
	var numThreads: Int = 4
	
	// How many threads can the tools use?
	// This will be passed to BWA and GATK tools.
	@Argument(
		doc="Integer number of jobs to scatter across nodes/cpus.",
		shortName="nsc",
		fullName="number_of_scatters",
		required=false
	)
	var numScatters: Int = 1
	
	// Path to bwa executable.
	// Should we pad /bwa on the end of it's missing?
	@Argument(
		doc="System path to bwa command.",
		shortName="bwa",
		fullName="bwa_path",
		required=false
	)
	var bwaCmd: File = _
	
	// Path to qualimap executable.
	@Argument(
		doc="System path to qualimap command.",
		shortName="qm",
		fullName="qualimap_path",
		required=false
	)
	var qualimapCmd: File = _
	
	@Argument(
		doc="Allow different capture platforms to merge.",
		shortName="xp",
		fullName="allow_cross_platform",
		required=false
	)
	val crossPlatform: Boolean = false
	
	@Argument(
		doc="HaplotypeCaller include mitochondria with the given ploidy. If less than 2 then effectively not specified as mitochondria will be excluded.",
		shortName="mp",
		fullName="mitochondira_ploidy",
		required=false
	)
	var mitoPloidy: Int = 0
	
	@Argument(
		doc="File Structure. 'single' to operate in individual's folder. Anything else will create subfolders for dna, library and runs",
		shortName="s",
		fullName="file_folder_structure",
		required=false
	)
	var fileStructure: String = "folder"
	
	val DEBUG: Boolean = false
	
	// Capture platforms enumeration
	// Collected list of capture platforms linked to their .bed file.
	var capturePath:	File				= null
	var capPlatforms:	Map[String,File]	= Map()	// PLATFORM -> FILEPATH.bed
	var capData:		Map[String,File]	= Map()	// PLATFORM -> FILEPATH.bed
	//var capBounds:		Map[String,String]	= Map()	// PLATFORM -> x,x,...
	
	// Environment
	//val tmpDir:String = System.getenv("HOME") + "/.queueTmp/"
	
	//////////////////
	// Base Configs //
	//////////////////
	
	//val settingsFile:	File = "settings.txt"
	
	// Generally known Indel sites to examine.
	var mills:	File = null
	var hapmap:	File = null
	var dbsnp:	File = null
	
	// Localized list of known sites to examine.
	var knownIndels:	Seq[File] = Nil
	
	// Process Files
	
	// Stage 1: Primary Alignment, sort and clean bam
	val paBWAAligned:	File = "paAligned"	//"all"
	val paSortedSam:	File = "paSorted"	//"allsort"
	val paCleanedSam:	File = "paCleaned"	//"allsortclean"
	
	// Stage 2: Per-sample fingerprinting for Identity validation.
	val idcHaplotyped:	File = "idcHaplotyped"
	val idcSelected:	File = "idcSelected"
	val idcPrinted:		File = "idcFingerprint"
	val idcGenotyped:	File = "idcGenotyped"
	val idcMatched:		File = "idcMatched"
	val idcMatchFail:	File = "idcMatch_FAILED"
	
	// Stage 3: Merge runs to library.
	val m2lTargetInts:	File = "mrTarget"
	val m2lIndelRealign:File = "mrRealigned"
	val m2lMarkDups:	File = "mrMarkDuped"
	val m2lQualiMapped:	File = "mrQualiMap"
	
	// Stage 4: Merge libraries to individual
	val m2iQualiMapped:	File = "mlQualiMap"
	val m2iTargetInts:	File = "mlTarget"
	val m2iIndelRealign:File = "mlRealigned"
	val m2iBQSR1:		File = "BQSR_FirstPass"
	val m2iBQSR2:		File = "BQSR_SecondPass"
	val m2iACoverants:	File = "AnalyzeCovariates"
	val m2iPrintReads:	File = "PrintReads"
	
	// Stage 5: Gender valiation
	val gvMatched:		File = "gvMatched"
	val gvDepth:		File = "gvDepthOfCoverage"
	
	// Stage 6: Fingerprinting
	val fpHaplotyped:	File = "fpHaplotyped"
	val fpSelected:		File = "fpSelected"
	val fpPrinted:		File = "FingerPrint"
	
	// Stage 7: Variant calling
	val vcAutosomeOut:	File = "vcAutosome"
	val vcPar1XOut:		File = "vcParOneX"
	val vcTrueXOut:		File = "vcTrueX"
	val vcPar2XOut:		File = "vcParTwoX"
	val vcFullYOut:		File = "vcFullY"
	val vcGLOut:		File = "vcGL"
	val vcMitoOut:		File = "vcMitochondria"
	val vcCatVar:		File = "vcCatVarred"
	
	// Shell scripts holders
	var fingerPrintCmd:		File = null
	var sampleMatchCmd:		File = null
	var bridgeLibraryCmd:	File = null
	var genderMatchCmd:		File = null
	
	// File Breaker
	// "/" for folder sub-division
	// "_" for underscore (single directory) sub-division
	var fileBreak: File = null	// "/"	// " + fileBreak + "
	
	// Grabbed once on first *.fastq.gz -> *_trimmed.fastq.gz
	// Re-parsed for additional files.
	// Saves re-walking over the WHOLE trimmed folder... which is big and slow!
	var trimmedFastQBlob: List[File] = Nil
	
	// Common variants for sample comparison
	var commonVariants: File = null
	
	var par1XRegion: String = null
	var trueXRegion: String = null 
	var par2XRegion: String = null
	var trueYRegion: String = null
	
	// Define the FastQdescriptions headers.
	// This is for reference and debugging purposes only.
	var fqHead: Seq[String] = Nil
	/*Seq(
		"Individual ID",		// 0
		"DNA Number",			// 1
		"Library Number",		// 2
		"Run Number",			// 3
		"Sequence Platform",	// 4
		"Filename1",			// 5
		"Filename2",			// 6
		"Directory Location",	// 7
		"Project Name",			// 8
		"Exome chip",			// 9
		"Target Coverage",		// 10
		"Date sequenced",		// 11
		"Family ID",			// 12
		"Individual ID",		// 13
		"Father ID",			// 14
		"Mother ID",			// 15
		"Gender",				// 16
		"Affection Status",		// 17
		"Any Other Comments",	// 18
		"Status"				// 19
	)*/
	
	// Frequently used header entry positions.
	var fqID:	Int = 0	// 0
	var fqDNA:	Int = 0	// 1
	var fqLib:	Int = 0	// 2
	var fqRun:	Int = 0	// 3
	var fqCap:	Int = 0	// 9
	var fqGen:	Int = 0	// 16
	var fqChr:	Int = 0	// 16
	
	// Where it all begins //
	
	def script = {
		PrintLog("----------\nFastQ to Genomic VCF\n----------\n", "AlignReport.txt")
	
		// If it exists, read settings file
		if(
			settingsFile.exists &&
			settingsFile.canRead
		) {
			grepConfigs(settingsFile)
			
			PrintLog("INFORMA\tSettings imported:" + 
				"\n\tAlignment request file: " + alignmentsFile +
				"\n\tSample description file: " + sampleDescriptionFile +
				"\n\tReference sequence file: " + referenceFile + 
				"\n\tBase output folder path: " + basePath + 
				"\n\tNumber of threads: " + numThreads + 
				"\n\tNumber of scatters: " + numScatters +
				"\n\tPath to BWA command: " + bwaCmd +
				"\n\tPath to Qualimap command: " + qualimapCmd +
				"\n\tCross platform merging is " + (if (!crossPlatform) "NOT ") + "permitted." +
				"\n\tMitochondria is " + (if (mitoPloidy > 2) "included at a ploidy of " + mitoPloidy else "NOT included") + "." +
				"\n\tFile structur is: " + fileStructure, "AlignReport.txt")
			
			if (knownIndels.size != knownIndels.distinct.size)
				PrintLog("WARNING\tKnown indels has duplicate entries!", "AlignReport.txt")
			
			grepPlatforms(capturePath)
			grepSource(capturePath)
			//grepPlatformBoundries(capturePath)
			
			// Check for missing sex chromosome tolerance boundries.
			for (platforms <- capPlatforms.keys)
				//if (!capBoundries.contains(platforms))
				if (!capData.contains(platforms))
					PrintLog("WARNING\tCapture platform " + platforms + " does not have gender bounds defined.", "AlignReport.txt")
			
			PrintLog("----------\n", "AlignReport.txt")
			
			// Set file generation mode.
			// Single: all files in single folder.
			// Anything else separates files by folder structure.
			// Already make a folder for the individual.
			fileBreak = if (fileStructure.toLowerCase == "single") "_" else "/"
			
			// Get alignment requests from request file.
			// Match alignment requests against entried in FastQ Description file. FQDF
			// Split samples by individual.
			grepIndividuals(grepFastQLines(grepSamples(alignmentsFile)))
		} else {
			PrintLog("\nFAILURE\tSettings file " + settingsFile + " does not exist!\n", "AlignReport.txt")
		}
	}
	
	/**
	 * Read core settings from file.
	 *
	 * Takes settings file location.
	 * Parses entries to set default paths, files and modes.
	 *
	**/
	
	def grepConfigs(configFile: File) {
		val setupLines = getFileLines(configFile)
		
		// Parse lines in settings file.
		for (line <- setupLines) {
			if (
				line.length > 0 &&		// Line has content
				!line.startsWith("#")	// Line doesn't start with a # for comments.
			) {
				// Splits based on any clump of white-space
				val pairs	= line.replaceAll("\"","").split("\\s+").toList
				
				val head:	String	= pairs(0).toLowerCase
				val param:	String	= pairs(1)
				
				println("INFORMA\tSetting " + head + ": " + param)
				
				if (head == "par1xregion")
					par1XRegion			= param
				else if (head == "truexregion")
					trueXRegion			= param
				else if (head == "par2xregion")
					par2XRegion			= param
				else if (head == "trueyregion")
					trueYRegion			= param
				else if (
					head == "basepath" &&
					basePath == null
				) basePath				= param
				else if (new File(param).exists) {
					// Display incoming values.
					if (
						head == "reference" &&
						referenceFile == null
					) referenceFile			= param
					else if (
						head == "sampledesc" &&
						sampleDescriptionFile == null
					) sampleDescriptionFile	= param
					else if (
						head == "common" &&
						commonVariants == null
					) commonVariants		= param
					else if (
						head == "capturepath" &&
						capturePath == null
					) capturePath			= param
					else if (
						head == "bwacmd" &&
						bwaCmd == null
					) bwaCmd				= param
					else if (
						head == "qualimapcmd" &&
						qualimapCmd == null
					) qualimapCmd			= param
					else if (head == "samplematchcmd")
						sampleMatchCmd		= param
					else if (head == "bridgelibcmd")
						bridgeLibraryCmd	= param
					else if (head == "fingerprintcmd")
						fingerPrintCmd		= param
					else if (head == "gendermatchcmd")
						genderMatchCmd		= param
					else if (head == "known")
						knownIndels		  :+= param
					else if (head == "mills")
						mills				= param
					else if (head == "hapmap")
						hapmap				= param
					else if (head == "dbsnp")
						dbsnp				= param
				} else {
					PrintLog("WARNING\tUnmapped settings line: " + line + ". If it's a file, does it exist?", "AlignReport.txt")
				}
			}
		}
		knownIndels = knownIndels.distinct
	}
	
	/**
	 * Parses the alignment request file for entries.
	 *
	 * Takes alignment request file.
	 * Returns sequence of alignment requests.
	 */
	
	def grepSamples(alignFile: File): Seq[String] = {
		// Ensure file exists.
		if (
			!alignFile.exists ||
			!alignFile.canRead
		) {
			PrintLog("\nFAILURE\tAlignment request file \"" + alignFile + "\" does not exist or cannot be read.\n", "AlignReport.txt")
			return Nil
		}
		
		// Read lines from alignments request file
		// Convert lines to a list
		// Sorted the list
		val alignList = getFileLines(alignFile).sorted
		
		//PrintLog("INFORMA\tAlignment file:", "AlignReport.txt")
		//for(requests <- alignList) PrintLog("\t" + requests, "AlignReport.txt")
		
		// add warning about duplication of alignment request.
		if (alignList.size != alignList.distinct.size) {
			// Print out any duplicates found.
			var misAlignment: String = ""
			for (alignment <- alignList) {
				var alignMatch: Int = 0
				for(alignAgain <- alignList)
					if (alignment.startsWith(alignAgain))
						alignMatch = alignMatch + 1
				
				if (alignMatch > 1)
					misAlignment = misAlignment +
						(if (misAlignment.length > 0) ", ") +
						alignment.replace("\\s+"," ")
			}
			PrintLog("\nFAILURE\tSamples list in \"" + alignFile + "\" contains duplicated entries for " + misAlignment + "! Please review sample list and fix!\n", "AlignReport.txt")
			
			return Nil
		}
		
		// Build list of samples.
		var samples: Seq[String] = Nil
		
		// Cycle through all entries in alignments list.
		for (line <- alignList) {
			// Split input entry by any continuous regions of whitespace.
			// Add Entires that where found in the fastqdescription file.
			if (
				!line.isEmpty &&
				!line.startsWith("#")
			) {
				val stitchLine = line.split("\\s+").toList
				samples :+= stitchLine(0) + '\t' +
							stitchLine(1) + '\t' +
							stitchLine(2) + '\t' +
							stitchLine(3)
			}
		}
		
		return samples
	}
	
	/** 
	 * Parses FastQ Description File for matching alignment requests
	 *
	 * Takes sequence of alignment requests
	 * Returns Sequence of lines from FQDF that matched alignment requests.
	 */
	
	def grepFastQLines(newAlign:Seq[String]): Seq[String] = {
		// Read source as fastq description file
		var alignList: Seq[String] = newAlign
		
		if (
			!sampleDescriptionFile.exists ||
			!sampleDescriptionFile.canRead
		) {
			PrintLog("\nFAILURE\tSample descriptions file \"" + sampleDescriptionFile + "\" does not exist or cannot be read! Got permissions?\n", "AlignReport.txt")
			
			return Nil
		}
		
		val descriptions: List[String] = getFileLines(sampleDescriptionFile)
		
		// validate for duplicate 1st 4 blocks in fqdf for error checking.
		// kill it if duplication found for current sample list
		// warn if dupilciate found that does not affect the requested runs
		
		var lineMatches: Seq[String] = Nil
		
		// Cycle through lines in fastq description file
		for (line <- descriptions) {
			val lineBlocks: Seq[String] = line.split("\t").toList
				
			if (lineBlocks.size > 4) {
				// Get header names for smart matching later.
				if (
					line.contains("Individual") &&
					line.contains("DNA") &&
					line.contains("Library") &&
					line.contains("Run")
				) {
					fqHead	= line.split("\\t").toList
					fqID	= fqHead.indexOf("Individual ID")
					fqDNA	= fqHead.indexOf("DNA Number")
					fqLib	= fqHead.indexOf("Library Number")
					fqRun	= fqHead.indexOf("Run Number")
					fqCap	= fqHead.indexOf("Exome chip")
					fqGen	= fqHead.indexOf("Gender")
					fqChr	= fqHead.indexOf("Sex Chromosomes")
				}
				
				// Cycle through alignment requests.
				for (alignRequest <- newAlign) {
					val alignBlocks: Seq[String] = alignRequest.split("\\s").toList
					
					// Scan 1st 4 positions from fqdf.
					// Can't use .startsWith since "a b c d" matches "a b c de"
					if (
						alignBlocks(0) == lineBlocks(0) &&
						alignBlocks(1) == lineBlocks(1) &&
						alignBlocks(2) == lineBlocks(2) &&
						alignBlocks(3) == lineBlocks(3)
					) {
						// this line matched a request. Add to the list.
						lineMatches :+= line
						
						PrintLog("INFORMA\tAlignement request " + alignRequest.replaceAll("\t",",") + " matched!", "AlignReport.txt")
						
						// Check if matches lines without current line is more than 1.
						if (lineMatches.filter(x => x == line).size > 1) {
							// More than 1 matched line in list. This means multiple entires in FQDF!
							PrintLog("\nFAILURE\tMultiple matches for " + alignRequest.replace("\t"," ") + " found in \"" + sampleDescriptionFile + "\"!\n\tPlease review before re-running pipeline!\n", "AlignReport.txt")
							
							return Nil
						}
						
						alignList = alignList.filterNot(x => x == alignRequest)
					}
				}
			} // else it's a blank line or a comment
		}
		
		if (alignList.nonEmpty) {
			// We have more requests than matched lines!
			// This mean there are unmatched request lines.
			PrintLog("\nFAILURE\t" + alignList.size + " unmatched samples:", "AlignReport.txt")
			for (request <- alignList) PrintLog("\n\t" + request.replace("\t"," "), "AlignReport.txt")
			PrintLog("", "AlignReport.txt")
			
			return Nil
		}
		
		// Returns sorted by individual, DNA, Library, Run, Platform, etc...
		return lineMatches
	}
	
	/**
	 * Locate a trimmed version of the given fastq file.
	 *
	 * Take path and file name and search for trimmed versions of a file.
	 * Returns the the trimmed or untrimmed file, in that order.
	 */
	
	def grepTrimmedFastQFile(path: File, file: File): File = {
		val rePath: File = path.split("/Originals")(0) + "/Trimmed"
		
		val previousRun: File = swapSuffix(alignmentsFile, ".txt", ".trimmed.list")
		
		if (previousRun.exists) {
			trimmedFastQBlob = getFileLines(previousRun)
		} else if (trimmedFastQBlob.size == 0) {
			// Trimmed trimmedFastQBlob not filled so grab all files in trimmed folder.
			PrintLog("INFORMA\tNo trimmed fastq file blob. Walking " + rePath, "AlignReport.txt")
			
			// This SUCKS, but there's really no other way.
			// It's a relative drop in the bucket.
			// Filter for *.fastq.gz
			// Filter for *_trimmed*
			// Filter out *unpaired*
			trimmedFastQBlob = grepFileTree(new File(rePath)).filter(_.getName.endsWith("fastq.gz")).filter(_.getName.contains("_trimmed")).filterNot(_.getName.contains("unpaired")).toList
			
			// Dump list of matched files to list.
			printToFile(previousRun) { p =>
				trimmedFastQBlob.foreach(p.println)
			}
		}
		
		// Locate my trimmed file from trimmed directory blob.
		
		// *.fastq.gz -> *._trimmed.fastq.gz -> *._trimmed_paired.fastq.gz
		val myBlob: List[File] = trimmedFastQBlob.filter(_.getName.startsWith(file.split(".fastq")(0)))
		
		if (myBlob.size == 1) {
			PrintLog("INFORMA\tFound trimmed file " + myBlob(0).getName, "AlignReport.txt")
			return myBlob(0)
		} else if (myBlob.size > 1) {
			PrintLog("\nFAILURE\tToo many matches for " + file.getName + " found in " + rePath + ":", "AlignReport.txt")
			for (files <- myBlob) println("\t" + files.getName)
			PrintLog("", "AlignReport.txt")
			
			return null
		}
		
		PrintLog("\nFAILURE\tNo matches for " + file.getName + " found in " + rePath + ".\n", "AlignReport.txt")
		return null
	}
	
	/**
	 * Gather all runs for an individual and pass them to primary alignment.
	 */
	
	def grepIndividuals(samples: Seq[String]) {
		// Replicate sample list for per-individual pruning.
		var sampleList: Seq[String] = samples
		
		// Cycle through all samples.
		for (sample <- sampleList) {
			// Split sample up to blocks.
			val sampleBlock: Seq[String] = sample.split("\t").toList
			
			// List of matched individuals for this cycle.
			var samplesGrepped: Seq[String] = Nil
			
			// Compare current cycle to list of samples.
			for (compare <- sampleList) {
				val compareBlock: Seq[String] = compare.split("\t").toList
				
				if (compareBlock(fqID) == sampleBlock(fqID)) {
					// Individual match!
					
					// Add compared sample to current individual.
					samplesGrepped :+= compare
					
					// Prune compare entry from list.
					sampleList = sampleList.filterNot(x => x == compare)
				}
			}
			
			// Check that more than 0 matches as can be 0 from pruned list iterating over pruned entries.
			if (samplesGrepped.size > 0) {
				// More than 0 matches found.
				val greppedBlock:	Seq[String]	= samplesGrepped(0).split("\t")
				val individual:		String		= greppedBlock(fqID)
				val greppedChip:	String		= greppedBlock(fqCap)
				val greppedGender:	String		= greppedBlock(fqGen)
				val greppedChromes:	String		= greppedBlock(fqChr)
				val greppedSample:	String		= greppedBlock(fqID) + "_" +
												  greppedBlock(fqDNA) + "_" +
												  greppedBlock(fqLib) + "_" +
												  greppedBlock(fqRun)
				var passChecks:		Boolean		= true
				
				PrintLog("----------\nMerge report for " + individual + "\n----------\n", basePath + "/" + greppedChip + "/" + individual + "/MergeReport.txt")
				
				PrintLog(
					"INFORMA\tID " + individual +
					" has " + samplesGrepped.size +
					" samples.", basePath + "/" + greppedChip + "/" + individual + "/MergeReport.txt"
				) // ** Output to individial/MergeReport.txt **
				
				// Verify consistency for gender and platform matching.
				// If only 1 sample, most wont matter.
				for (consistencyCheck <- samplesGrepped) {
					val checkBlock:		Seq[String]	= consistencyCheck.split("\t").toList
					val checkChip:		String = checkBlock(fqCap)
					val checkGender:	String = checkBlock(fqGen)
					val checkChromes:	String = checkBlock(fqChr)
					val checkSample:	String = checkBlock(fqID) + "_" +
												 checkBlock(fqDNA) + "_" +
												 checkBlock(fqLib) + "_" +
												 checkBlock(fqRun)
					
					// Gender consistency from fqdf.
					if (
						checkGender == "" ||
						checkGender == "0"
					)	// Gender not specified.
						PrintLog("WARNING\tGender is not specified for sample " + checkSample + "\n\tIf all samples for this individual are unknown then gender will be automatically mapped later.", basePath + "/" + checkChip + "/" + individual + "/MergeReport.txt")

					if (checkGender != greppedGender)						// Gender varies between two samples for an individual.
						PrintLog("WARNING\tGender for sample " + checkSample + ": " + checkGender + " does not match sample " + greppedSample + ": " + greppedGender + "\n\tPlease review sample consistency check output " + "to make sure samples are from the same individual.", basePath + "/" + checkChip + "/" + individual + "/MergeReport.txt")
					
					if (checkChromes != greppedChromes)
						PrintLog("WARNING\tSex Chromosomes for sample " + checkSample + ": " + checkChromes + " does not match sample " + greppedSample + ": " + greppedChromes + "\n\tPlease review sample consistency check output " + "to make sure samples are from the same individual.", basePath + "/" + checkChip + "/" + individual + "/MergeReport.txt")
					
					if (checkChip != greppedChip) {
						PrintLog(
							"WARNING\tSample " + checkSample +
							" Exome chip does not match " + greppedSample
						, basePath + "/" + checkChip + "/" + individual + "/MergeReport.txt")
						if (crossPlatform) {
							PrintLog("WARNING\tHowever, cross-platform merging is allowed.", basePath + "/" + checkChip + "/" + individual + "/MergeReport.txt")
						} else {
							PrintLog("\nFAILURE\tCross-platform merging is NOT allowed!\n", "AlignReport.txt")
							passChecks = false
						}
					}
					
					if (
						checkChip == "" ||
						checkChip == "unknown"
					) {
						PrintLog("WARNING\tExome chip is unknown.\n\tPlease review qualimap output and adjust.\n\tThis pipeline will stop after PrintReads", basePath + "/Unknown/" + individual + "/MergeReport.txt")
					}
				}
				
				// Did all samples pass gender and platform matching?
				if (passChecks) {
					// Yes! Split samples by Library for given individual
					PrimaryAlignSortCleanValidate(samplesGrepped)
					grepLibraries(samplesGrepped)
					MergeLibrariesToIndividual(samplesGrepped)
				} else {
					// No. Warn about this individual.
					PrintLog("\nFAILURE\tExcluding Individual " + individual + " from pipeline!\n\tPlease verify samples and sample descriptions file for errors.\n", "AlignReport.txt")
					
					println("\tPress Control + C to cancel run.\n\tContinuing with remaining individuals in 60 seconds.\n")
					
					Thread.sleep(if (DEBUG) 10000 else 60000)
				}
			}
		}
	}
	
	/**
	 * Splits samples up by library for a given individual.
	 *
	 * Takes a list of samples for an individual.
	 * Passes runs to grepRuns to process runs.
	 * Calls MergeRunsToLibrary on current list of libraries.
	 *
	 * This code should look kind of familiar :(
	 */
	
	def grepLibraries(samples:Seq[String]) {
		// replicate sample list for pruning.
		var sampleList: Seq[String] = samples
		
		// Cycle through libraries for this individual
		for (sample <- sampleList) {
			// Break up sample line blocks
			val sampleBlock = sample.split("\t").toList
			
			// Matched library container
			var samplesGrepped: Seq[String] = Nil
			
			for (compare <- sampleList) {
				// Break up compare line blocks
				val compareBlock = compare.split("\t").toList
				
				if (
					compareBlock(fqID) == sampleBlock(fqID) &&	// Re-check individual match. You neven know who'll edit what.
					compareBlock(fqLib) == sampleBlock(fqLib)
				) {
					// Individual and Library match!
					
					// Add found sample to list.
					samplesGrepped :+= compare
					
					// Prune compared sample from list.
					sampleList = sampleList.filterNot(x => x == compare)
				}
			}
			
			// Check for matched libraries still in list.
			if (samplesGrepped.size > 0) {
				// At least 1 library matching
				val greppedBlock:	Seq[String]	= samplesGrepped(0).split("\t")
				val individual:		String		= greppedBlock(fqID)
				val library:		String		= greppedBlock(fqLib)
				val samplesFound:	Int			= samplesGrepped.size
				val greppedChip:	String		= greppedBlock(fqCap)
				val greppedSample:	String		= greppedBlock(fqID) + "_" +
												  greppedBlock(fqDNA) + "_" +
												  greppedBlock(fqLib) + "_" +
												  greppedBlock(fqRun)
				var passChecks:		Boolean		= true
				
				PrintLog(
					"INFORMA\tID " +
					individual +
					" Library " +
					library +
					" has " + samplesFound + " sample" + (if (samplesFound != 1) "s") + " to merge."
				, basePath + "/" + greppedChip + "/" + individual + "/MergeReport.txt")
				
				// Verify consistency for gender and platform matching.
				// If only 1 sample, most wont matter.
				for (consistencyCheck <- samplesGrepped) {
					val checkBlock:		Seq[String] = consistencyCheck.split("\t").toList
					val checkChip:		String = checkBlock(fqCap)
					val checkSample:	String = checkBlock(fqID) + "_" +
												 checkBlock(fqDNA) + "_" +
												 checkBlock(fqLib) + "_" +
												 checkBlock(fqRun)
					
					if (checkChip != greppedChip) {
						PrintLog(
							"\nFAILURE\tSample " + checkSample +
							" Exome chip: " + checkChip + " does not match " + greppedSample + ": " + greppedChip + " within a library!\n"
						, "AlignReport.txt")
						passChecks = false
					} else if (
						checkChip == "" ||
						checkChip == "unknown"
					) {
						PrintLog("WARNING\tExome chip is unknown.\n\tPlease review qualimap output and adjust.\n\tThis pipeline will stop after PrintReads", basePath + "/" + checkChip + "/" + individual + "/MergeReport.txt")
					}
				}
				
				// Did all samples pass gender and platform matching?
				if (passChecks) {
					// Yes! Split samples by Library for given individual
					MergeRunsToLibrary(samplesGrepped)
				} else {
					// No. Warn about this individual.
					PrintLog("\nFAILURE\tExcluding Individual " + individual + " from pipeline!" +
						"\n\tPlease verify samples and sample descriptions file for errors.\n", "AlignReport.txt")
						
					println("FAILURE\tPress Control + C to cancel run." +
						"\n\tContinuing with remaining individuals in 60 seconds.\n")
					
					Thread.sleep(if (DEBUG) 10000 else 60000)
				}
			}
		}
	}
	
	/**
	 * Runs BWA Mem, SortSam and CleanSam (3n) on a given sample.
	 *
	 * Takes a line from the Sample Descriptions File.
	 *
	 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 * !!! Will there ever be mulitple runs that are from different cap platforms? !!!
	 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 */
	
	def PrimaryAlignSortCleanValidate(samples:Seq[String]){
		// Get the base sample to compare stuff to, if there's more than 1.
		val zeroBlock:		Seq[String]	= samples(0).split("\t").toList
		val individual:		String		= zeroBlock(fqID)
		val dnaZero:		String		= zeroBlock(fqDNA)
		val libraryZero:	String		= zeroBlock(fqLib)
		val runZero:		String		= zeroBlock(fqRun)
		val genderZero:		String		= zeroBlock(fqGen)
		val chromeZero:		String		= zeroBlock(fqChr)
		val capPlatZero:	String		= zeroBlock(fqCap)
		
		// Cycle through all samples.
		for (sample <- samples) {
			val bwa					= new BurrowsWheelerAligner
//			val bwa					= new BurrowsWheelerAlignerPicardSortSam
			val sortSam				= new SortSam				with Picard_Arguments
			val cleanSam			= new CleanSam				with Picard_Arguments
		
			val haplotypeCaller		= new HaplotypeCaller		with GATK_Arguments with ISR_Intersect
			val selectVariants		= new SelectVariants		with GATK_Arguments
			val fingerPrinter		= new FingerPrinter
			val genotypeConcordance	= new GenotypeConcordance	with GATK_Arguments
			
			val sampleMatcher		= new SampleMatcher
			
			val sampleBlock = sample.split("\t").toList
			
			// Get base filename.
			val fastqpath:	File	= sampleBlock(fqHead.indexOf("Directory Location"))
			
			// Seek out trimmed versions of these files.
			val filename1: File = 
				grepTrimmedFastQFile(
					fastqpath,
					sampleBlock(fqHead.indexOf("Filename1"))
				)
				
			val filename2: File =
				grepTrimmedFastQFile(
					fastqpath,
					sampleBlock(fqHead.indexOf("Filename2"))
				)
			
			if (
				filename1 != null &&
				filename2 != null
			) {
				
				// Read a block from fastq file.
				// Capture everything until 1st new-line char
				// Capture everything after 1st @ char
				// Capture @ delimited blocks.
				val headerBlocks:	Seq[String]	= zcat(filename1).split("\n")(0).split("@")(1).split(":").toList
				
				// May be required when machines names start with FC? Discuss!
				val numblocks:		Int		= headerBlocks.size
				
				//println("headerblocks " + numblocks)
				//headerBlocks.foreach(println)
				
				val dnaCur:			String	= sampleBlock(fqDNA)
				val libraryCur:		String	= sampleBlock(fqLib)
				val runCur:			String	= sampleBlock(fqRun)
				val genderCur:		String	= sampleBlock(fqGen)
				val chromeCur:		String	= sampleBlock(fqChr)
				val capPlatCur:		String	= sampleBlock(fqCap)
				val platform:		String	= sampleBlock(fqHead.indexOf("Sequence Platform"))
				
				val instrument:		String	= headerBlocks(0).replaceAll(" ","_")
				val instrumentrun:	String	= headerBlocks(1)
				val flowcell:		String	= headerBlocks(2)
				val lane:			String	= headerBlocks(3)
				val index:			String	= headerBlocks(9)
				
				// Construct output paths.
				val runPath:		String	= basePath +
										"/" + capPlatCur +
										"/" + individual + 
									"/DNA_" + dnaCur +
					 fileBreak + "Library_" + libraryCur +
						 fileBreak + "Run_" + runCur +
							 fileBreak
							 
				val runPathZero:	String	= basePath +
										"/" + capPlatZero +
										"/" + individual + 
									"/DNA_" + dnaZero +
					 fileBreak + "Library_" + libraryZero +
						 fileBreak + "Run_" + runZero +
							 fileBreak
				
				// Create output path if it doesn't already exist.
				if (runPath.contains("/"))
					new File(runPath.substring(0,runPath.lastIndexOf('/'))).mkdirs
				
				  ///////////////////////////////////
				 // Begin Primary Alignment phase //
				///////////////////////////////////
				
				// Burrows Wheeler Aligner FastQ x2 to SAM
				bwa.reference_sequence	= referenceFile.replace(".fasta","")
				bwa.operation			= "mem -M"
				bwa.input				= Seq(filename1,filename2)
				bwa.output				= runPath + paBWAAligned + ".sam"
				bwa.readGroups			= "@RG" + "'\\t'" + 
					"ID:" + instrument + "_" + instrumentrun + "_" +
							flowcell + "_" + lane + "_" + index + "'\\t'" + 
					"PL:" + platform + "'\\t'" + 
					"PU:" + flowcell + "." + lane + "'\\t'" + 
					"LB:" + individual + "_" + dnaCur + "_" + libraryCur + "'\\t'" + 
					"SM:" + individual
				
				bwa.nt					= numThreads
				//bwa.scatterCount		= 4
				
				bwa.jobName				= individual + "_PA1_BWA"
				bwa.analysisName		= individual + "_" +
										  dnaCur + "_" +
										  libraryCur + "_" +
										  runCur
				bwa.isIntermediate		= true
				
				// Sort SAM to BAM file
				sortSam.input			= Seq(bwa.output)
				sortSam.output			= runPath + paSortedSam + ".bam"
				sortSam.maxRecordsInRam	= 2000000
				
				sortSam.jobName			= individual + "_PA2_PSS"
				sortSam.analysisName	= bwa.analysisName
				
				// Clean BAM file & generate index
//				cleanSam.input			= Seq(bwa.output)
				cleanSam.input			= Seq(sortSam.output)
				cleanSam.output			= runPath + paCleanedSam + ".bam"
				
				cleanSam.jobName		= individual + "_PA3_PCS"
				cleanSam.analysisName	= bwa.analysisName
				
				  /////////////////////////////////
				 // End Primary Alignment phase //
				/////////////////////////////////
				
				  ////////////////////////////////
				 // Begin Identity Check phase //
				////////////////////////////////
				
				// Haplotype Called on specifid intervals. BAM to VCF
				haplotypeCaller.input_file		= Seq(cleanSam.output)
				haplotypeCaller.out				= runPath + idcHaplotyped + ".vcf"
				
				haplotypeCaller.genotyping_mode	=
					org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
					
				haplotypeCaller.alleles			= commonVariants
				haplotypeCaller.intervals		=
					Seq(
						commonVariants,
						capPlatforms(capPlatCur.toLowerCase)
					)
				
				haplotypeCaller.nct				= numThreads
				
				haplotypeCaller.analysisName	= bwa.analysisName
				haplotypeCaller.jobName			= individual + "_IDC1_HTC"
				
				fingerPrinter.input			= haplotypeCaller.out
				fingerPrinter.output		= runPath + idcSelected + ".vcf"
				
				fingerPrinter.analysisName	= bwa.analysisName
				fingerPrinter.jobName		= individual + "_IDC2_GFP"
				fingerPrinter.isIntermediate= true
				
				// Select individual variants from de-printed.
				selectVariants.variant		= fingerPrinter.output
				selectVariants.select		=
					Seq("vc.getGenotype('Identity').getGQ() > 10")
					
				selectVariants.out			= runPath + idcPrinted + ".vcf.gz"
				
				selectVariants.nt			= numThreads
				
				selectVariants.analysisName	= bwa.analysisName
				selectVariants.jobName		= individual + "_IDC3_SV"
				
				// Compare genotypes.
				genotypeConcordance.comp			= runPathZero + idcPrinted + ".vcf.gz"
				genotypeConcordance.eval			= selectVariants.out
				genotypeConcordance.out				= runPath + idcGenotyped + ".grp"
				genotypeConcordance.ignoreFilters	= true
				genotypeConcordance.moltenize		= true
				
				genotypeConcordance.analysisName	= bwa.analysisName
				genotypeConcordance.jobName			= individual + "_IDC4_GTC"
				
				// Extract genotype value from concordance output.
				sampleMatcher.ident				= basePath + "/" + capPlatCur + "/" + individual
				sampleMatcher.bamFile			= cleanSam.output
				sampleMatcher.concordanceFile	= genotypeConcordance.out
				sampleMatcher.truthFile			= genotypeConcordance.comp
				sampleMatcher.linkTo			= runPath + idcMatched + ".bam"
				
				if (
					genderCur != genderZero ||		// Gender match
					chromeCur != chromeZero ||		// Sex Chromosomes match
					(
						!crossPlatform &&			// Cross platform merging allowed
						capPlatCur != capPlatZero	// Capture platforms match
					)
				) {
					if (genderCur != genderZero)
						PrintLog("WARNING\tGender does not match in Sample Description file!", basePath + "/" + capPlatCur + "/" + individual + "/MergeReport.txt")
					
					if (chromeCur != chromeZero)
						PrintLog("WARNING\tSex Chromosomes do not match in Sample Description file!", basePath + "/" + capPlatCur + "/" + individual + "/MergeReport.txt")
					
					if (capPlatCur != capPlatZero)
						PrintLog("WARNING\tCapture platforms do not match in Sample Description file!", basePath + "/" + capPlatCur + "/" + individual + "/MergeReport.txt")
					
					PrintLog("FAILURE\tThis individual's pipeline will stop at gender validation.\n\t", basePath + "/" + capPlatCur + "/" + individual + "/MergeReport.txt")
					
					println("\tPlease press Control + C to cancel now.\n\tRemaining pipeline will continue in 60 seconds.")
					
					Thread sleep (if (DEBUG) 10000 else 60000)
				}
				
				sampleMatcher.analysisName		= bwa.analysisName
				sampleMatcher.jobName			= individual + "_IDC5_SM"
				sampleMatcher.isIntermediate	= true	// Keep output for testing purposes.
				
				  //////////////////////////////
				 // End Identity Check phase //
				//////////////////////////////
				
				add(
					bwa,
					sortSam
				)
				
				if (samples.size > 1) {
					// Multiple samples. compare them!
					add(
						cleanSam,
						haplotypeCaller,
						selectVariants,
						fingerPrinter,
						genotypeConcordance,
						sampleMatcher
					)
				} else {
					// Only 1 run for this individual so it matches by default.
					PrintLog("INFORMA\tOnly 1 run for " + individual + " so skipped self-matching results.", basePath + "/" + capPlatCur + "/" + individual + "/MergeReport.txt")
					cleanSam.output = runPath + idcMatched + ".bam"
					cleanSam.isIntermediate = true	// Keep output for testing purposes.
					add(cleanSam)
				}
				
			} else {
				PrintLog("\nFAILURE\tSuitable fastQ file missing!", "AlignReport.txt")
				PrintLog("FAILURE\t1: \"" + filename1 + "\"", "AlignReport.txt")
				PrintLog("FAILURE\t2: \"" + filename2 + "\"\n", "AlignReport.txt")
			}
		}
	}
	
	
	/**
	 * Merge all runs into the current library.
	 *
	 * Takes a sequence of lines from FastQ Description file.
	 * Returns true if Runs can be merged to a Library.
	 * Returns false if Runs cannot be merged, due to gender discrepency and
	 *  whatever other issues can be detected.o
	 */
	
	def MergeRunsToLibrary(samples:Seq[String]) {
		var individual:	String		= null
		var libPath:	String		= null
		var jobTitle:	String		= null
		var capture:	Seq[String]	= Nil
		var platform:	String		= null
		
		val realignerTargetCreator	= new RealignerTargetCreator	with GATK_Arguments
		val indelRealigner			= new IndelRealigner			with GATK_Arguments
		val markDuplicates			= new MarkDuplicates			with Picard_Arguments
		
		// Cycle through all passed samples to build list of BAM files for realignment.
		// Checks for gender matching across all bam files generated.
		for (sample <- samples) {
			// Get the blocks... yet again.
			val sampleLine: List[String] = sample.split("\t").toList
			
			// Useful info!
			val dna:		String	= sampleLine(fqDNA)
			val library:	String	= sampleLine(fqLib)
			val run:		String	= sampleLine(fqRun)
			
			individual	= sampleLine(fqID)
			platform	= sampleLine(fqCap)
			libPath		= "/DNA_" + dna + fileBreak + "Library_" + library + fileBreak
			jobTitle	= individual + "_" + dna + "_" + library
			
			if (
				(
					platform == "" ||
					platform.toLowerCase == "unknown"
				) &&
				capture.size < capPlatforms.keys.size
			) {
				// Sample has unknown capture platform
				// This hasn't already triggered for this library.
				PrintLog("\nFAILURE\tUnknown capture platform for this library." +
					"\n\tQualiMap will run for all platforms on this library." +
					"\n\tThis library will not be handed to next phase." +
					"\n\tThis library output will be appended with _unknown" +
					"\n\tThis individual's pipeline will fail." +
					"\n\tPlease review qualimap output for best-guess capture platform.\n" +
					"\n\tOnce you have determined capture platform remove _unknown from .bam and .metric file" +
					"\n\tThen re-run pipeline for this individual to pick up where you left off.\n", basePath + "/Unknown/" + individual + "/MergeReport.txt")
				
				println("FAILURE\tPress Control + C to cancel pipeline now" +
					"\n\tPipeline will continue in 60 seconds.\n")
				
				Thread.sleep(if (DEBUG) 10000 else 60000)
				
				capture = capPlatforms.keys.toList
			}
			
			// What overall benefit is there to including the BAM files here?
			// w/BAM inputs the outputs of this and the Indel Realignment files take
			//  longer to generate and are larger in size. If there in a benefit,
			//  does it require larger data sets?
			//  RTC Time and file size are almost doubled.
			//  IReal Time is more than doubled but file size if only slightly larger.
			
			realignerTargetCreator.input_file :+= basePath + "/" + platform + "/" + individual + "/" + libPath + "Run_" + run + fileBreak + idcMatched + ".bam"
		}
		
		// Ensure no duplicates bam files.
		capture = capture.distinct
		realignerTargetCreator.input_file	= realignerTargetCreator.input_file.distinct
		
		  /////////////////////////////////////////
		 // Generate indel realignment targets. //
		/////////////////////////////////////////
		
		realignerTargetCreator.known		= knownIndels
		println(capture.size + " " + platform)
		realignerTargetCreator.out			= 
			if (capture.size > 1)
				basePath + "/Unknown/" + individual + "/" + libPath + m2lTargetInts + ".intervals"
			else
				basePath + "/" + platform + "/" + individual + "/" + libPath + m2lTargetInts + ".intervals"
		
		// Parallelization
		realignerTargetCreator.nt			= numThreads
		realignerTargetCreator.scatterCount	= parallelThreadsAndScatters
		
		realignerTargetCreator.jobName		= individual + "_M2L_RTC"
		realignerTargetCreator.analysisName	= jobTitle
		
		  /////////////////////////////////////
		 // Library level Indel Realignment //
		/////////////////////////////////////
		
		indelRealigner.targetIntervals	= realignerTargetCreator.out
		indelRealigner.input_file		= realignerTargetCreator.input_file
		indelRealigner.out				= 
			if (capture.size > 1)
				basePath + "/Unknown/" + individual + "/" + libPath + m2lIndelRealign + ".bam"
			else
				basePath + "/" + platform + "/" + individual + "/" + libPath + m2lIndelRealign + ".bam"
		
		// Parallelization.
		indelRealigner.scatterCount		= parallelScatterOnly
		
		indelRealigner.jobName			= individual + "_M2L_IR"
		indelRealigner.analysisName		= jobTitle
		
		  /////////////////////
		 // Mark Duplicates //
		/////////////////////
		
		//markDuplicates.input			= Seq(indelRealigner.out)
		markDuplicates.input			= realignerTargetCreator.input_file
		markDuplicates.output			=
			if (capture.size > 1)
				basePath + "/Unknown/" + individual + "/" + libPath + m2lMarkDups + "_unknown" + ".bam"
			else
				basePath + "/" + platform + "/" + individual + "/" + libPath + m2lMarkDups + ".bam"
		
		markDuplicates.outputIndex		= swapSuffix(markDuplicates.output,".bam",".bai")
		markDuplicates.metrics			= swapSuffix(markDuplicates.output,".bam",".metric")
		
		markDuplicates.assumeSorted		= true
		
		markDuplicates.jobName			= individual + "_M2L_MD"
		markDuplicates.analysisName		= jobTitle
		markDuplicates.isIntermediate	= true	// Keep output for testing purposes.
		
		add(
			markDuplicates
		)
		//	realignerTargetCreator,
		//	indelRealigner,
		
		  //////////////////////////////////////////////////////
		 // Run qualimap if the capture platform was unknown //
		//////////////////////////////////////////////////////
		
		/*if (capture.size > 1) {
			for (
				platform <- capture
				if platform != "genomic"	// Skip the genomic QM.
			) {
				add(
					QualiMapper(
						bamFile = markDuplicates.output,
						outDir = basePath + "/Unknown/" + individual + "/" + libPath + m2lQualiMapped + "_" + platform,
						outFile = basePath + "/Unknown/" + individual + "/" + libPath + m2lQualiMapped + "_" + platform,
						numThreads = numThreads,
						platform = platform.toLowerCase,
						jobName = individual + "_L_QM",
						analysisName = jobTitle
					)
				)
			}
		}*/
	}
	
	/**
	 * Merge Library runs to individual.
	 *
	 */
	
	def MergeLibrariesToIndividual(allSamples:Seq[String]) {
		var libPath:		String		= null
		var jobTitle:		String		= null
		var indiPath:		String		= null
		var gender:			String		= null
		var sexChromosome:	String		= null
		var individual:		String		= null
		var capture:		Seq[String]	= Nil
		
		// Merging Libraries to individual
		val realignerTargetCreator	= new RealignerTargetCreator	with GATK_Arguments
		val indelRealigner			= new IndelRealigner			with GATK_Arguments
		
		// Finalise Individual BAM
		val baseRecalibrationFirst	= new BaseRecalibrator			with GATK_Arguments
		val baseRecalibrationSecond	= new BaseRecalibrator			with GATK_Arguments
		val analyzeCovariates		= new AnalyzeCovariates			with GATK_Arguments
		val printReads				= new PrintReads				with GATK_Arguments
		
		// Gender validation.
		val depthOfACoverage		= new DepthOfCoverage			with GATK_Arguments with ISR_Intersect with DoC_Arguments
		val depthOfXCoverage		= new DepthOfCoverage			with GATK_Arguments with ISR_Intersect with DoC_Arguments
		val depthOfYCoverage		= new DepthOfCoverage			with GATK_Arguments with ISR_Intersect with DoC_Arguments
		val genderMatch				= new GenderMatch
		
		// Finger printer
		val haplotypeCaller			= new HaplotypeCaller			with GATK_Arguments with ISR_Intersect
		val selectVariants			= new SelectVariants			with GATK_Arguments
		val fingerPrinter			= new FingerPrinter
		
		// Variant caller
		// autosomeHtC is done separately.
		val glHtC		= new HaplotypeCaller	with GATK_Arguments with VariantIndex_Arguments with VCHtC_Arguments
		val par1XHtC	= new HaplotypeCaller	with GATK_Arguments with VariantIndex_Arguments with VCHtC_Arguments
		val trueXHtC	= new HaplotypeCaller	with GATK_Arguments with VariantIndex_Arguments with VCHtC_Arguments
		val par2XHtC	= new HaplotypeCaller	with GATK_Arguments with VariantIndex_Arguments with VCHtC_Arguments
		val fullYHtC	= new HaplotypeCaller	with GATK_Arguments with VariantIndex_Arguments with VCHtC_Arguments
		val mitoHtC		= new HaplotypeCaller	with GATK_Arguments with VariantIndex_Arguments with VCHtC_Arguments
		
		val catVar			= new CatVariants //org.broadinstitute.gatk.tools.CatVariants
		catVar.reference	= glHtC.reference_sequence
		//catVar.assumeSorted	= false	// Only if autosomal is not per-chromosome.
		
		  /////////////////////////
		 // Job Paralellization //
		/////////////////////////
		
		realignerTargetCreator.nt			= numThreads
		realignerTargetCreator.scatterCount	= parallelThreadsAndScatters()
		indelRealigner.scatterCount			= parallelScatterOnly()
		
		baseRecalibrationFirst.nct			= numThreads
		baseRecalibrationFirst.scatterCount	= realignerTargetCreator.scatterCount
		
		baseRecalibrationSecond.nct			= numThreads
		baseRecalibrationSecond.scatterCount= realignerTargetCreator.scatterCount
		
		analyzeCovariates.scatterCount		= indelRealigner.scatterCount
		
		printReads.nct						= numThreads
		printReads.scatterCount				= realignerTargetCreator.scatterCount
		
		//haplotypeCaller.nct					= numThreads
		selectVariants.nt					= numThreads
		  //////////////////////
		 // Output retention //
		//////////////////////
		
		analyzeCovariates.isIntermediate	= false	// ALWAYS keep this .PDF
		printReads.isIntermediate			= false	// ALWAYS keep this .BAM
		selectVariants.isIntermediate		= false	// ALWAYS keep this .VCF.GZ
		
		
		  /////////////////////////////////////////////////////
		 // Build list of libraries to merge if applicable. //
		/////////////////////////////////////////////////////
		
		for (sample <- allSamples) {
			val sampleLine: List[String] = sample.split("\t").toList
			
			println(s"${sampleLine}")
			
			individual				= sampleLine(fqID)
			val dna:		String	= sampleLine(fqDNA)
			val library:	String	= sampleLine(fqLib)
			var platform:	String	= sampleLine(fqCap)
			printf("DEBUG\tGender [%s]:", gender)
			gender					=
				if (gender != "") {			// Already defined as undefined
					if (gender == null) {	// Initial value.
						sampleLine(fqGen)	// Get stated value
					} else if (gender != sampleLine(fqGen)) {	// Doesn't match previous instance.
						""	// Define as undefined.
					} else {	// Matches previous instance and isn't null.
						gender	// Retain value.
					}
				} else {	// Is Undefined.
					gender	// Retain value.
				}
			printf("[%s]\n", gender)
									
			printf("DEBUG\tChrome [%s]:", sexChromosome)
			sexChromosome			=
				if (	// Sometimes the sex chromosome column is missing. Lets not freak out about that...
					fqChr > -1 &&
					fqChr < sampleLine.size
				) {
					if (sexChromosome != "") {				// Already defined as undefined
						if (sexChromosome == null) {		// Initial value.
							sampleLine(fqChr).toLowerCase	// Get stated value
						} else if (sexChromosome != sampleLine(fqChr).toLowerCase) {	// Doesn't match previous instance.
							""	// Define as undefined.
						} else {			// Matches previous instance and isn't null.
							sexChromosome	// Retain value.
						}
					} else {	// is Underfined.
						sexChromosome	// Retain value.
					}
				} else {	// Unable to retrieve sex chromosome column in sample description file.
					""		// Define as undefined.
				}
			printf("[%s]\n", sexChromosome)
			
			indiPath	= basePath + "/" + platform + "/" + individual + "/"
			libPath		= indiPath + "DNA_" + dna + fileBreak + "Library_" + library + fileBreak
			
			capture :+= platform
			
			realignerTargetCreator.input_file :+= libPath + m2lMarkDups + ".bam"
		}
		
		  /////////////////////
		 // Job identifiers //
		/////////////////////
		
		realignerTargetCreator.analysisName	= individual
		indelRealigner.analysisName			= individual
		
		baseRecalibrationFirst.analysisName	= individual
		baseRecalibrationSecond.analysisName= individual
		analyzeCovariates.analysisName		= individual
		printReads.analysisName				= individual
		
		haplotypeCaller.analysisName		= individual
		fingerPrinter.analysisName			= individual
		selectVariants.analysisName			= individual
		
		depthOfACoverage.analysisName		= individual
		depthOfXCoverage.analysisName		= individual
		depthOfYCoverage.analysisName		= individual
		genderMatch.analysisName			= individual
		
		par1XHtC.analysisName				= individual
		trueXHtC.analysisName				= individual
		par2XHtC.analysisName				= individual
		fullYHtC.analysisName				= individual
		glHtC.analysisName					= individual
		mitoHtC.analysisName				= individual
		catVar.analysisName					= individual
		
		realignerTargetCreator.jobName		= individual + "_M2I_RTC"
		indelRealigner.jobName				= individual + "_M2I_IR"
		
		baseRecalibrationFirst.jobName		= individual + "_BR1"
		baseRecalibrationSecond.jobName		= individual + "_BR2"
		analyzeCovariates.jobName			= individual + "_ACV"
		printReads.jobName					= individual + "_PR"
		
		haplotypeCaller.jobName				= individual + "_FP_HTC"
		fingerPrinter.jobName				= individual + "_FP"
		selectVariants.jobName				= individual + "_FP_SV"
		
		depthOfACoverage.jobName			= individual + "_IGV_DoAC"
		depthOfXCoverage.jobName			= individual + "_IGV_DoXC"
		depthOfYCoverage.jobName			= individual + "_IGV_DoYC"
		genderMatch.jobName					= individual + "_IGV_Match"
		
		par1XHtC.jobName					= individual + "_VC_PAR1X"
		trueXHtC.jobName					= individual + "_VC_TrueX"
		par2XHtC.jobName					= individual + "_VC_PAR2X"
		fullYHtC.jobName					= individual + "_VC_TrueY"
		glHtC.jobName						= individual + "_VC_GL"
		mitoHtC.jobName						= individual + "_VC_Mito"
		catVar.jobName						= individual + "_VC_CatVar"
		
		// Prune duplicate entries. This occurs due to incoming list being sample based instead of library based.
		capture = capture.distinct
		realignerTargetCreator.input_file = realignerTargetCreator.input_file.distinct
		
		// Only realign if more than one library to merge.
		if (realignerTargetCreator.input_file.size > 1) {
			// more than one library to merge.
			
			PrintLog("INFORMA\tID " + individual + " has " + realignerTargetCreator.input_file.size + " libraries to merge!", indiPath + "MergeReport.txt")
			realignerTargetCreator.input_file.foreach(println)
			
			  /////////////////////////////////////////////////
			 // Individual level Realignment Target Creator //
			/////////////////////////////////////////////////
			
			realignerTargetCreator.known		= knownIndels
			realignerTargetCreator.out			= indiPath + m2iTargetInts + ".intervals"
			
			  ////////////////////////////////////////////
			 // Individual level Realignment and Merge //
			////////////////////////////////////////////
			
			indelRealigner.input_file		= realignerTargetCreator.input_file
			indelRealigner.targetIntervals	= realignerTargetCreator.out
			indelRealigner.out				= indiPath + m2iIndelRealign + ".bam"
			
			add(
				realignerTargetCreator,
				indelRealigner
			)
			
			// No repeating markDuplicates at individual level?
			// What if there are duplicates across libraries?
		
		} else {
			// there was only 1 library to merge so lets skip re-aligning it!
			
			PrintLog("INFORMA\tID " + individual + " only has 1 library! No merging required.", indiPath + "MergeReport.txt")
			
			indelRealigner.out = libPath + m2lMarkDups + ".bam"
		}
		
		  //////////////////////////////////////
		 // Base Quality Score Recalibration //
		//////////////////////////////////////
		
		baseRecalibrationFirst.knownSites	= Seq(mills,dbsnp)
		baseRecalibrationFirst.input_file	= Seq(indelRealigner.out)
		baseRecalibrationFirst.out			= indiPath + m2iBQSR1 + ".table"
		
		  //////////////////////
		 // BQSR Second Pass //
		//////////////////////
		
		baseRecalibrationSecond.knownSites	= baseRecalibrationFirst.knownSites
		baseRecalibrationSecond.input_file	= Seq(indelRealigner.out)
		baseRecalibrationSecond.BQSR		= baseRecalibrationFirst.out
		baseRecalibrationSecond.out			= indiPath + m2iBQSR2 + ".table"
		
		  ////////////////////////
		 // Analyze Covariates //
		////////////////////////
		
		analyzeCovariates.before		= baseRecalibrationFirst.out
		analyzeCovariates.after			= baseRecalibrationSecond.out
		analyzeCovariates.plots			= indiPath + m2iACoverants + ".pdf"
		
		  /////////////////
		 // Print Reads //
		/////////////////
		
		printReads.input_file	= Seq(indelRealigner.out)
		printReads.BQSR			= baseRecalibrationFirst.out
		printReads.out			= indiPath + individual + "." + m2iPrintReads + ".bam"
		
		add(
			baseRecalibrationFirst,
			printReads
		)
		
		//	baseRecalibrationSecond,
		//	analyzeCovariates,
		
		val platform:	String	= capture(0).toLowerCase
		
		if (capture.size == 1) {
			// Only 1 capture platform.
			// Since we didn't run qualimap per-library, run it now.
			
			add(
				QualiMapper(
					bamFile		= printReads.out,
					outDir		= indiPath + "/" + m2iQualiMapped + "_" + platform,
					outFile		= indiPath + "/" + m2iQualiMapped + "_" + platform,
					numThreads	= numThreads,
					platform	= platform.toLowerCase,
					jobName		= individual + "_QM",
					analysisName= individual
				)
			)
		} else {
			PrintLog("INFORMA\tQuality Mapping not done at individual level for " + individual + " as there were multiple platforms involved." +
				"\n\tPlease review library level quality maps for reference as required.", indiPath + "MergeReport.txt")
		}
		
		  ////////////////////////
		 // Global Fingerprint //
		////////////////////////
		
		// Haplotype Called on specifid intervals. BAM to VCF
		haplotypeCaller.input_file	= Seq(printReads.out)
		haplotypeCaller.out			= indiPath + "/" + idcHaplotyped + ".vcf"
		
		haplotypeCaller.genotyping_mode	= org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
		
		haplotypeCaller.alleles		= commonVariants
		haplotypeCaller.dbsnp		= dbsnp
		haplotypeCaller.intervals	= Seq(commonVariants,capPlatforms(platform))
		
		// Replace identity with "Individual"
		fingerPrinter.input			= haplotypeCaller.out
		fingerPrinter.output		= indiPath + "/" + idcSelected + ".vcf"
		
		// Select individual variants from vcf.
		selectVariants.variant		= fingerPrinter.output
		//selectVariants.variant		= haplotypeCaller.out
		selectVariants.select		=
			Seq("vc.getGenotype('Identity').getGQ() > 30")
//			Seq("vc.getGenotype('" + individual + "').getGQ() > 30")
			
		selectVariants.out			= indiPath + "/" + individual + "_" + idcPrinted + ".vcf.gz"
		
		add(
			haplotypeCaller,
			fingerPrinter,
			selectVariants
		)
		
		  ///////////////////////
		 // Gender Validation //
		///////////////////////
		
		// Run gender check.
		val depthOfCoverageOutput	= indiPath + "/" + gvDepth
		
		depthOfACoverage.input_file	= Seq(printReads.out)
		depthOfXCoverage.input_file	= depthOfACoverage.input_file
		depthOfYCoverage.input_file	= depthOfACoverage.input_file
		
		// Because Genomic repeat sequences in Y chromosome causes issues.
		depthOfACoverage.intervals	= if (platform == "Genomic") {
										Seq(capPlatforms("AV5"))
									} else {
										Seq(capPlatforms(platform))
									}
		depthOfXCoverage.intervals	= depthOfACoverage.intervals
		depthOfYCoverage.intervals	= depthOfACoverage.intervals
		
		depthOfACoverage.intervalsString= Seq("1")
		depthOfXCoverage.intervalsString= Seq(trueXRegion)
		depthOfYCoverage.intervalsString= Seq(trueYRegion)
		
		// It actually outputs .sample_summary and .sample_statistic files.
		// This means this step technically produces no output and will be
		// repeated every time you run this individual.
		depthOfACoverage.out	= depthOfCoverageOutput + "-A"
		depthOfXCoverage.out	= depthOfCoverageOutput + "-X"
		depthOfYCoverage.out	= depthOfCoverageOutput + "-Y"
		
		genderMatch.ident		= indiPath
		genderMatch.gender		= gender
		genderMatch.sChro		= sexChromosome.toUpperCase
		//genderMatch.threshold	= capBounds(platform).replaceAll("\t",",")
		genderMatch.threshold	= capData(platform)
		genderMatch.xChro		= depthOfXCoverage.out
		genderMatch.yChro		= depthOfYCoverage.out
		genderMatch.aChro		= depthOfACoverage.out
		genderMatch.bamFile		= printReads.out
		genderMatch.bamLink		= indiPath + "/" + gvMatched + ".bam"
		
		// If gender doesn't match output will not be generated?
		// Does this preserve this function's input files even if they are
		// marked as intermediate by creator?
		
		add(
			depthOfXCoverage,
			depthOfYCoverage,
			depthOfACoverage,
			genderMatch
		)
		
		  /////////////////////
		 // Variant Calling //
		/////////////////////
		
		// If we get here without a sex chromosome
		// in the sample description file, lets assume.
		if (
			sexChromosome == "" ||
			sexChromosome == "unknown" ||
			sexChromosome == "undefined"
		) {
			if (gender == "1") sexChromosome = "xy"
			else if (gender == "2") sexChromosome = "xx"
		}
		
		val sexChromosomeXCount:Int		= sexChromosome.count('x'==)
		val sexChromosomeYCount:Int		= sexChromosome.count('y'==)
		
		var tmpFolder: File = qSettings.tempDirectory
		
		for (Chromosome <- 1 to 22) {
			val autosomeHtC				= new HaplotypeCaller with GATK_Arguments with VariantIndex_Arguments with VCHtC_Arguments
			autosomeHtC.analysisName	= individual
			autosomeHtC.jobName			= individual + "_VC_Chr" + Chromosome
			autosomeHtC.input_file		= Seq(genderMatch.bamLink)
			autosomeHtC.intervalsString	= Seq(""+Chromosome)	// Cuz they're stringy.
			autosomeHtC.sample_ploidy	= 2
			
			autosomeHtC.out				=
				if (autosomeHtC.scatterCount > 1)
					tmpFolder + "/" + individual + "/" + vcAutosomeOut + "_" + Chromosome + ".g.vcf"
				else
					indiPath + "/" + vcAutosomeOut + "_" + Chromosome + ".g.vcf"
					
			catVar.variant :+= autosomeHtC.out
			
			add(autosomeHtC)
		}
		
		// Inputs
		par1XHtC.input_file		= Seq(genderMatch.bamLink)
		trueXHtC.input_file		= Seq(genderMatch.bamLink)
		par2XHtC.input_file		= Seq(genderMatch.bamLink)
		fullYHtC.input_file		= Seq(genderMatch.bamLink)
		mitoHtC.input_file		= Seq(genderMatch.bamLink)
		glHtC.input_file		= Seq(genderMatch.bamLink)
		
		// Intervals
		par1XHtC.intervalsString	= Seq(par1XRegion)
		trueXHtC.intervalsString	= Seq(trueXRegion)
		par2XHtC.intervalsString	= Seq(par2XRegion)
		fullYHtC.intervalsString	= Seq("Y")
		mitoHtC.intervalsString		= Seq("MT")
		glHtC.excludeIntervalsString= Seq("1","2","3","4","5","6","7","8","9",
			"10","11","12","13","14","15","16","17","18","19","20","21","22",
			"X","Y","MT")
		
		// Ploidy
		glHtC.sample_ploidy			= 2
		
		// Calulate ploidy based on actual sex chromosome counts.
		par1XHtC.sample_ploidy		= if (
										sexChromosomeXCount == 2 &&
										sexChromosomeYCount == 2
									)
										2
									else
										(sexChromosomeXCount + sexChromosomeYCount)
		
		trueXHtC.sample_ploidy	= sexChromosomeXCount
		par2XHtC.sample_ploidy	= par1XHtC.sample_ploidy
		
		// Always include a Y read. This will insert no-read data for females.
		// 99.8+% of time this will be 1. But there's a 1:500+ chance for XYY
		// or something.
		fullYHtC.sample_ploidy	= if (sexChromosomeYCount > 1)
									sexChromosomeYCount
								else
									1
		// Only implimented when mitoPloidy is > 2
		mitoHtC.sample_ploidy	= mitoPloidy
		
		// Outputs
		// If autosomal is scattered, output final file to local tmp folder.
		// This uses CatVariants witch drags down when catting networked files.
		// If single file, output to main directory.
		par1XHtC.out			= indiPath + "/" + vcPar1XOut + ".g.vcf"
		trueXHtC.out			= indiPath + "/" + vcTrueXOut + ".g.vcf"
		par2XHtC.out			= indiPath + "/" + vcPar2XOut + ".g.vcf"
		fullYHtC.out			= indiPath + "/" + vcFullYOut + ".g.vcf"
		mitoHtC.out				= indiPath + "/" + vcMitoOut + ".g.vcf"
		glHtC.out				= indiPath + "/" + vcGLOut + ".g.vcf"
		
		// Output catVar to local storage
		// Network storage can start to grind.
		catVar.out				= tmpFolder + "/" + individual + "/" + individual + ".g.vcf.gz"
		
		// Build catvar input
		catVar.variant :+= par1XHtC.out
		catVar.variant :+= trueXHtC.out
		catVar.variant :+= par2XHtC.out
		catVar.variant :+= fullYHtC.out
		catVar.variant :+= glHtC.out
		
		add(
			par1XHtC,
			trueXHtC,
			par2XHtC,
			fullYHtC,
			glHtC
		)
		
		if (mitoPloidy > 1) {
			PrintLog("INFORMA\tIncluding mitochondria in this pipeline with a ploidy of " + mitoPloidy, indiPath + "MergeReport.txt")
			add(mitoHtC)
			catVar.variant :+= mitoHtC.out
		} else {
			PrintLog("INFORMA\tExcluding mitochondria from this pipeline!", indiPath + "MergeReport.txt")
		}
		
		add(catVar)
		
		// Move temp file to main storage.
		add(new CommandLineFunction {
			@Input(doc="incoming file",required=true)
			var input: File		= catVar.out
			
			@Output(doc="outgoing file",required=true)
			var outFile: File	= indiPath + "/" + individual + ".copy"
			
			this.isIntermediate	= false
			this.jobName		= individual + "_VC_MoveTo"
			this.analysisName	= individual
			
			def commandLine =
				"/home/clinical_genetics/bin/q-shell/movefiles.sh '" + input +
				"' '" + outFile + "'"
		})
	}
	
	/**
	 * Base arguments.
	 */
	
	trait GATK_Arguments extends CommandLineGATK {
		this.reference_sequence	= referenceFile
		this.nCoresRequest		= numThreads
		this.isIntermediate		= true
	}
	
	trait ISR_Intersect extends CommandLineGATK {
		this.isr = org.broadinstitute.gatk.utils.interval.IntervalSetRule.INTERSECTION
	}
	
	trait DoC_Arguments extends DepthOfCoverage {
		this.omitDepthOutputAtEachBase	= true
		this.omitLocusTable				= true
		this.omitIntervals				= true
	}
	
	trait VariantIndex_Arguments extends CommandLineGATK {
		this.variant_index_type		= org.broadinstitute.gatk.utils.variant.GATKVCFIndexType.LINEAR
		this.variant_index_parameter= Some(128000)
	}
	
	trait VCHtC_Arguments extends HaplotypeCaller {
		//this.dbsnp				= dbsnp
		this.emitRefConfidence	= org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF
		this.pairHMM			= org.broadinstitute.gatk.utils.pairhmm.PairHMM.HMM_IMPLEMENTATION.LOGLESS_CACHING
		//this.nct				= numThreads	// Do not set -nct for HaplotypeCaller, it does nothing but eat up memory.
	}
	
	trait Picard_Arguments extends PicardBamFunction {
		this.nCoresRequest		= numThreads
		this.isIntermediate		= true
	}
	
	/**
	 * Output useful number of scatters to use.
	 */ 
	 
	def parallelScatterOnly(): Int = if (numScatters > 1 && numScatters > numThreads) return numScatters else return numThreads
	
	def parallelThreadsAndScatters(): Int = if (numScatters > 1) return numScatters else return (numThreads / 2)
	
	/**
	 * Return File descriptor with oldExt exchanged for newExt
	 */
	
	def swapSuffix(file: File, oldExt:String, newExt: String): File = new File(file.getAbsolutePath.stripSuffix(oldExt) + newExt)
	
	/**
	 * Returns all lines from a file as a list of strings.
	 */
	
	def getFileLines(file: File): List[String] = scala.io.Source.fromFile(file).getLines.toList
	
	/**
	 * Build list of capture platforms
	 */
	
	def grepPlatforms(path: File) {
		val fileBlob = grepFileTree(new File(path)).filter(_.getName.endsWith(".bed"))
		for (file <- fileBlob)
			capPlatforms = 
				capPlatforms + 
				(file.getName().stripSuffix(".bed").toLowerCase -> file)
		
		println("INFORMA\tCapture platforms:")
		for (platform <- capPlatforms)
			println("\t" + platform)
	}
	
	/**
	 * Build list of capture platforms
	 */
	
	def grepSource(path: File) {
		val fileBlob = grepFileTree(new File(path)).filter(_.getName.endsWith(".sh"))
		for (file <- fileBlob)
			capData = 
				capData + 
				(file.getName().stripSuffix(".sh").toLowerCase -> file)
		
		println("INFORMA\tCapture boundries:")
		for (platform <- capData)
			println("\t" + platform)
	}
	
	/**
	 * Build list of capture platform gender tolerance boundries.
	 */
	/*
	def grepPlatformBoundries(path: File) {
		val fileBlob = grepFileTree(new File(path)).filter(_.getName.endsWith(".tsv"))
		for (file <- fileBlob) {
			val fileLines = getFileLines(file)
			for (line <- fileLines)
				if (
					!line.isEmpty &&
					!line.startsWith("#")
				) capBounds =
					capBounds +
					(file.getName().stripSuffix(".tsv").toLowerCase -> line)
		}
	}
	*/
	/**
	 * Write to a file
	 */
	 
	def PrintLog(string: String, file: File = "MergeReport.txt", append: Boolean = true) {
		
		println(string)	// Print to StdOut
		
		if (file.contains("/"))
			new File(file.substring(0,file.lastIndexOf('/'))).mkdirs
		
		val outFile = scala.tools.nsc.io.File(file)
		
		if (append && outFile.canWrite)
			outFile.appendAll(string + "\n")	// Log to file.
		else
			outFile.writeAll(string + "\n")		// Log to file.
	}
	
	/**
	 * Dump array/list/etc to file
	 */
	
	def printToFile(f: java.io.File)(op: java.io.PrintWriter => Unit) {
		val p = new java.io.PrintWriter(f)
		try { op(p) } finally { p.close() }
	}
	
	/**
	 * Recursive directory tree search that takes a relatively small amount of time!
	 *
	 * Takes a directory path.
	 * Returns a file stream containing all files within the base folder.
	 */
	
	def grepFileTree(f: File): Stream[File] = f #:: (
		if (f.isDirectory) {
			f.listFiles().toStream.flatMap(grepFileTree) 
		} else {
			Stream.empty
		}
	)
	
	/**
	 * 'zcat' command line replacement.
	 *
	 * Takes a gzipped file.
	 * Returns string containing the first chunk of bytes.
	 */
	
	def zcat(source: File): String = {
		// Buffer to store bytes from file.
		val buf = new Array[Byte](1024)
		
		// Attempt to retrieve bytes from file.
		try {
			// Get buf length number of bytes from source file.
			val in = new java.util.zip.GZIPInputStream(new FileInputStream(source)).read(buf)
			
			// Convert stream of bytes to a string.
			val str = new String(buf, "UTF-8")
			
			return str
		} catch {	// ERRORS? Errors!
			case _:FileNotFoundException =>
				System.err.printf("File Not Found: %s", source)
			case _:SecurityException =>
				System.err.printf("Permission Denied: %s", source)
		}
		return null
	}
	
	/**
	 * BWA mem command line job.
	 *
	 * Takes a pair of fastq files
	 * Outputs an aligned SAM file
	 *
	 * Constructs unsorted aligned SAM file from trimmed paired fastq files.
	 * Does this primary alignment help the indel realignment that comes later?
	 * Does indel realignment have to re-sort things?
	 * Can we just pipe BWA MEM output directly to Picard's SortSam?
	 * Could we use samtools' SortSam instead?
	 */
	
	class BurrowsWheelerAligner extends CommandLineFunction {
		jobName			= "BurrowsWheelerAligner"
		analysisName	= jobName
		
		@Input(
			doc="Reference sequence to align to.",
			required = true
		)
		var reference_sequence: File = _
		
		@Input(
			doc="The input fastq files to run primary alignment on.",
			required = true
		)
		var input: Seq[File] = Nil
		
		@Output(
			doc="The sorted BAM or SAM output file.",
			required = true
		)
		var output: File = _
		
		@Argument(
			doc="Number of threads to allocate to BWA",
			required = false
		)
		var nt: Int = _
		
		@Argument(
			doc="Read Group output header information.",
			required = false
		)
		var readGroups: String = _
		
		@Argument(
			doc="Operating parameters, algorith, etc",
			required=true
		)
		var operation: String = _
		
		@Argument(
			doc="Log file output location other than default",
			required=false
		)
		var logout: File = _
		
		override def commandLine =
			bwaCmd +
			required(operation, escape=false) +
			optional("-t", nt, escape=false) +
			optional("-R", readGroups, escape=false) +
			required(reference_sequence, escape=false) +
			repeat(input) +
			required(">", output, escape=false) + 
			optional("2>", logout, escape=false)
	}
	
	/**
	 * Grab the Overall_Genotype_Concordance from GenotypeConcordance function.
	 *
	 * Reads the GenotypeConcordance out files for the match ratio
	 * If the match ratio and alleles matched count are sufficient, then samples
	 *  share identity and the the pipeline is connected.
	 */
	
	class SampleMatcher extends CommandLineFunction {
		jobName			= "SampleMatcher"
		analysisName	= jobName
		
		@Argument(
			doc="Individual's identity string",
			required=true
		)
		var ident: String = _
		
		@Input(
			doc="Incoming file to be linked on successful matching.",
			required=true
		)
		var bamFile: File = _
		
		@Input(
			doc="Genotype concordance file to extract matching from.",
			required=true
		)
		var concordanceFile: File = _
		
		@Input(
			doc="Genotype concordance comp/truth set file.",
			required=true
		)
		var truthFile: File = _
		
		@Output(
			doc="Outgoing hardlinked file to incoming file.",
			required=true
		)
		var linkTo: File = _
		
		@Argument(
			doc="Minimum number of alleles matched by genotype concordance.",
			required=false
		) 
		var minAlleles: Int = 500
		
		@Argument(
			doc="Minimum matching ratio below which the linkage fails.",
			required=false
		)
		var minRatio: Double = 0.85
		
		var baiFile: File = null
		var baiLink: File = null
		
		override def freezeFieldValues() {
			super.freezeFieldValues()
			if (baiFile == null) {
				baiFile = swapSuffix(bamFile,".bam",".bai")
			}
			if (baiLink == null && linkTo != null) {
				baiLink = swapSuffix(linkTo,".bam",".bai")
			}
		}
		
		override def commandLine =
			required(sampleMatchCmd) + 
			required("-d", ident) +
			optional("-a", minAlleles) +
			required("-b", bamFile) + 
			required("-i", baiFile) +
			required("-c", concordanceFile) +
			required("-l", linkTo) + 
			required("-o", baiLink) +
			optional("-r", minRatio) +
			required("-t", truthFile)
	}
	
	/**
	 * Replace Individual ID number in Select variant fingerprint file
	 * This will allow comparison of any 2 samples regardless of actual identity.
	 *
	 * takes source file name
	 * takes output file name
	 * takes generic id. default: FingerPrint
	 *
	 * This output should be retained forever.
	 */
	
	class FingerPrinter extends CommandLineFunction {
		jobName			= "FingerPrinter"
		analysisName	= jobName
		isIntermediate	= true
		
		@Input(
			doc="Fingerprint input .vcf file.",
			required=true
		)
		var input: File = _
		
		@Output(
			doc="Fingerprint output .vcf file.",
			required=true
		)
		var output: File = _
		
		@Argument(
			doc="Fingerprint identity.",
			required=false
		)
		var genericID: String = "Identity"
		
		override def commandLine =
			required(fingerPrintCmd) +
			required("-i", input) +
			required("-f", genericID) +
			required("-o", output)
	}
	
	/**
	 * Qualitymap wrapper since we seem to do this in several spots.
	 */ 
	
	def QualiMapper(
		bamFile: File,
		outDir: File,
		outFile: File = null,
		numThreads: Int,
		platform: String,
		jobName: String = null,
		analysisName: String = null
	): QualiMap = {
		val qualityMap						= new QualiMap
		
		qualityMap.bam						= bamFile
		qualityMap.outdir					= outDir
		qualityMap.outfile					= 
			if ( outFile == null ) platform + ".pdf"
			else outFile + ".pdf"
		
		qualityMap.headless					= true
		
		qualityMap.number_of_threads		= numThreads
		
		if (platform.toLowerCase == "genomic") {
			qualityMap.memory				= "24G"
			qualityMap.number_of_reads		= 500
		} else {
			qualityMap.feature_file			= capPlatforms(platform.toLowerCase)
			qualityMap.memory				= "4G"
			qualityMap.outside_stats		= true
		}
		
		qualityMap.genome_distribution		= "HUMAN"
		qualityMap.paint_chromosome_limits	= true
		
		if (jobName != null)
			qualityMap.jobName				= jobName
			
		if (analysisName != null)
			qualityMap.analysisName			= analysisName
		
		return qualityMap
	}
	
	/**
	 * Qualimap command line job
	 *
	 */
	 
	class QualiMap extends CommandLineFunction {
		jobName			= "QualiMap"
		analysisName	= jobName
		isIntermediate	= false
		
		@Input(
			doc="BAM file input",
			shortName="bam",
			fullName="bam",
			required=true
		)
		var bam: File = _
		
		@Argument(
			doc="Paint chromosome limits inside charts",
			shortName="c",
			fullName="paint-chromosome-limits",
			required=false
		)
		var paint_chromosome_limits: Boolean = false
		
		@Argument(
			doc="Species to compare with genome GC distribution. Possible values: HUMAN or MOUSE.",
			shortName="gd",
			fullName="genome-gc-distr",
			required=false
		)
		var genome_distribution: String = "HUMAN"
		
		@Input(
			doc="Feature file with regions of interest in GFF/GTF or BED format",
			shortName="gff",
			fullName="feature-file",
			required=false
		)
		var feature_file: File = _
		
		@Argument(
		doc="Minimum size for a homopolymer to be considered in indel analysis (default is 3)",
			shortName="hm",
			fullName="homopolymer",
			required=false
		)
		var homopolymer: Int = _
		
		@Argument(
			doc="Activate this option to collect statistics of overlapping paired-end reads",
			shortName="ip",
			fullName="collect-overlap-pairs",
			required=false
		)
		var collect_overlap_pairs: Boolean = false
		
		@Argument(
			doc="Number of reads analyzed in a chunk (default is 1000)",
			shortName="nr",
			fullName="nr",
			required=false
		)
		var number_of_reads: Int = _
		
		@Argument(
			doc="Number of threads (default is 8)",
			shortName="nt",
			fullName="nt",
			required=false
		)
		var number_of_threads: Int = _
		
		@Argument(
			doc="Number of windows (default is 400)",
			shortName="nw",
			fullName="nw",
			required=false
		)
		var number_of_windows: Int = _
		
		@Argument(
			doc="File to save per base non-zero coverage. WARN  large files are  expected for large genomes",
			shortName="oc",
			fullName="output-genome-coverage",
			required=false
		)
		var output_genome_coverage: File = _
		
		@Argument(
			doc="Report information for the regions outside those defined by feature-file  (ignored when -gff option is not set)",
			shortName="os",
			fullName="outside-stats",
			required=false
		)
		var outside_stats: Boolean = false
		
		@Argument(
			doc="Output folder for html.",
			shortName="outdir",
			fullName="outdir",
			required=false
		)
		var outdir: File = _
		
		@Output(
			doc="Output file for PDF.",
			shortName="outfile",
			fullName="outfile",
			required=false
		)
		var outfile: File = _
		
		@Argument(
			doc="Output format. html or pdf.",
			shortName="outformat",
			fullName="outformat",
			required=false
		)
		var outformat: String = "pdf"
		
		@Argument(
			doc="Sequencing library protocol: strand-specific-forward, strand-specific-reverse or non-strand-specific (default)",
			shortName="p",
			fullName="sequencing-protocol",
			required=false
		)
		var sequencing_protocol: String = _
		
		@Argument(
			doc="Activate this option to skip duplicated alignments from analysis. If the duplicates are not flagged in BAM file, then they will be detected by Qualimap.",
			shortName="sd",
			fullName="skip-duplicated",
			required = false
		)
		var skip_duplicated: Boolean = false
		
		@Argument(
			doc="Java VM's memory allocation.",
			shortName="m",
			fullName="memory",
			required = false
		)
		var memory: String = _
		
		@Argument(
			doc="Java VM run as headless.",
			shortName="hl",
			fullName="headless",
			required = false
		)
		var headless: Boolean = true
		
		var fileOut: String = null
		
		// Need to strip the folder structure cuz qualimap dev are funny!
		override def freezeFieldValues() {
			super.freezeFieldValues()
			fileOut = outfile.getName
		}
		
		override def commandLine = required(qualimapCmd) +
			required("bamqc") +
			required("-bam", bam) +
			required("-outfile", fileOut) + 
			conditional(paint_chromosome_limits, "-c") + 
			conditional(collect_overlap_pairs, "-ip") + 
			conditional(skip_duplicated, "-sd") +
			conditional(outside_stats, "-os") + 
			//conditional(headless, "--Djava.awt.headless=true") + 
			optional("--java-mem-size=", memory, spaceSeparated=false) +
			optional("-outformat", outformat) + 
			optional("-outdir", outdir) + 
			optional("-gff", feature_file) + 
			optional("-gd", genome_distribution) + 
			conditional(homopolymer > 0, "-hm " + homopolymer) + 
			conditional(number_of_reads > 0, "-nr " + number_of_reads) +
			optional("-nt", number_of_threads) +
			conditional(number_of_windows > 0, "-nw " + number_of_windows) +
			optional("-oc", output_genome_coverage) + 
			optional("-p", sequencing_protocol)
	}
	
	/**
	 * Matches coverage of sex chromosomes to determine chromosomal gender.
	 */
	
	class GenderMatch extends CommandLineFunction {
		jobName			= "GenderMatch"
		analysisName	= jobName
		isIntermediate	= true
		
		@Argument(
			doc="Individual's Identity",
			required=true
		)
		var ident: String = _
		
		@Argument(
			doc="Gender column from sample description.",
			required=true
		)
		var gender:	String = _
		
		@Argument(
			doc="Sex Chromosome column from sample description.",
			required=false
		)
		var sChro:	String = _
		
		@Argument(
			doc="File to pass for per-platform X/Y:A ratio bias.",
			required=true
		)
		var threshold:	File = _
		
		@Input(
			doc="Depth of Coverage from X chromosome.",
			required=true
		)
		var xChro:	File = _
		
		@Input(
			doc="Depth of Coverage from Y chromosome.",
			required=true
		)
		var yChro:	File = _
		
		@Input(
			doc="Depth of Coverage from Autosomal chromosome.",
			required=true
		)
		var aChro:	File = _
		
		@Input(
			doc="Incoming .BAM file to be linked to.",
			required=true
		)
		var bamFile:File = _
		
		@Input(
			doc="Incoming .BAM file to be linked to.",
			required=false
		)
		var baiFile:File = _
		
		@Output(
			doc="Outgoing .BAM file coming from the link.",
			required=true
		)
		var bamLink:File = _
		
		@Output(
			doc="Outgoing .BAM file coming from the link.",
			required=false
		)
		var baiLink:File = _
		
		override def freezeFieldValues() {
			super.freezeFieldValues()
			if (baiFile == null && bamFile != null)
				baiFile = swapSuffix(bamFile,".bam",".bai")
			
			if (baiLink == null && bamLink != null)
				baiLink = swapSuffix(bamLink,".bam",".bai")
		}
		
		override def commandLine =
			required(genderMatchCmd) + 
			required("-d", ident) + 
			required("-x", xChro) +
			required("-y", yChro) +
			required("-a", aChro) +
			required("-b", bamFile) +
			required("-i", baiFile) +
			required("-l", bamLink) +
			required("-o", baiLink) +
			required("-p", threshold) +
			conditional((
				gender == "1" ||
				gender == "2"
			),"-g " + gender, escape=false) + 
			conditional((
				sChro != "" &&
				sChro != null
			), "-s " + sChro, escape=false)
	}
	
	/*
	class DepthOfCoverage_File extends DepthOfCoverage {
		override def freezeFieldValues() {
			super.freezeFieldValues()
			output = new File(output.getAbsolutePath.stripSuffix(".sample_summary"))
		}
	}*/
	
	/**
	 * Rename files!
	 *
	 * Moved a file from input to output.
	 *
	 * Used exclusively to bridge the library level merge to the individual.
	 * Run when cross-platform merging is allowed or when only 1 platform per library.
	 * Not running this will explicitly fail further jobs that require the given output.
	 */
	
	class BridgeLibraryPipeline extends CommandLineFunction {
		jobName			= "BridgePipe"
		analysisName	= jobName
		
		@Input(
			doc="Incoming BAM file to be linked.",
			required=true
		)
		var input: File = _
		
		// This is technically required.
		// We allow this to be blank to break the pipeline under specific circumstances
		@Output(
			doc="Outgoing file to be linked.",
			required=false
		)
		var output: File = _
		
		var baiFile: File = null
		var baiLink: File = null
		
		override def freezeFieldValues() {
			super.freezeFieldValues()
			baiFile = swapSuffix(input,".bam",".bai")
			if (output != null) {
				baiLink = swapSuffix(output,".bam",".bai")
			}
		}
		
		override def commandLine =
			required(bridgeLibraryCmd) + 
			required("-b", input) + 
			required("-i", baiFile) +
			optional("-l", output) +
			optional("-o", baiLink)
	}
	
}

/**
 * Picard CleanSam command line job.
 *
 * Take aligned and sorted BAM file.
 * Outputs a cleaned BAM file.
 */
 
class CleanSam extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction with PicardBamFunction {
	jobName			= "CleanSam"
	analysisName	= jobName
	javaMainClass	= "picard.sam.CleanSam"
	
	@Input(
		doc="The input SAM or BAM files to clean.",
		required = true
	)
	var input: Seq[File] = Nil
	
	@Output(
		doc="The sorted BAM or SAM output file.",
		required = true
	)
	var output: File = _
	
	@Output(
		doc="The output bam index",
		required = false
	)
	var outputIndex: File = _
	
	sortOrder = null
	
	override def freezeFieldValues() {
		super.freezeFieldValues()
		if (outputIndex == null && output != null) {
			outputIndex = new File(output.getAbsolutePath.stripSuffix(".bam") + ".bai")
		}
	}
	
	override def inputBams = input
	override def outputBam = output
	this.createIndex = Some(true)
	override def commandLine = super.commandLine
}

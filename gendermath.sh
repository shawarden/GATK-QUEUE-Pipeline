#!/bin/bash

# This scrip will compare the coverage output to determine chromosomal gender.
# 
# This is used to validate genders and identities.
#
# @params
#  i  Input VCF file
#  f  Generic Identity (defaults to "FingerPrint"
#  o  Output VCF file
#

HARDLINKER="/home/clinical_genetics/bin/q-shell/hardlink.sh"

usage() {
cat << EOF
This scrip will parse a vcf file and replace the individual ID

 Once generated it can be used to compare any other samples from this platform.

usage: $0 options

Options:
  -g  base gender specified
  -s  hard coded gender
  -p  platform file
  -x  X chromosome coverage file
  -y  Y chromosome coverage file
  -a  Autosomal chromosome coverage file
  -b  Source .BAM file
  -i  Source .BAI file
  -l  Linked .BAM file if gender match passes
  -o  Linked .BAI file if gender match passes
EOF
}

GENDER=Unknown
SEXCHR=UNKNOWN
XCOVER=
YCOVER=
ACOVER=
SRCBAM=
SRCBAI=
LNKBAM=
LNKBAI=
IDENTITY=

while getopts "a:b:d:g:i:l:o:s:p:x:y:" OPTION
do
	FILE=
	case $OPTION in
		d)
			IDENTITY=${OPTARG}
			;;
		g)
			if [[ "${OPTARG}" == "1" ]]
			then
				GENDER="Male"
			elif [[ "${OPTARG}" == "2" ]]
			then
				GENDER="Female"
			else
				echo "INFORMA: Base gender invalid." | tee -a "${IDENTITY}/MergeReport.txt"
				GENDER="Unknown"
			fi
			;;
		s)
			if [[ -z ${OPTARG} ]]
			then
				echo "INFORMA: No sex chromosomes defined. Comparing to base Gender." | tee -a "${IDENTITY}/MergeReport.txt"
				SEXCHR="Undefined"
			else
				SEXCHR=${OPTARG}
			fi
			;;
		x)
			FILE="${OPTARG}.sample_summary"
			if [[ -e ${FILE} ]]
			then
				echo "" >> ${OPTARG}	# make a fake output file for DoC
				XCOVER=${FILE}
			else
				echo "FAILURE: X coverage file ${FILE} does not exist." | tee -a "${IDENTITY}/MergeReport.txt"
			fi
			;;
		y)
			FILE="${OPTARG}.sample_summary"
			if [[ -e ${FILE} ]]
			then
				echo "" >> ${OPTARG}	# make a fake output file for DoC
				YCOVER=${FILE}
			else
				echo "FAILURE: Y coverage file ${FILE} does not exist." | tee -a "${IDENTITY}/MergeReport.txt"
			fi
			;;
		a)
			FILE="${OPTARG}.sample_summary"
			if [[ -e ${FILE} ]]
			then
				echo "" >> ${OPTARG}	# make a fake output file for DoC
				ACOVER=${FILE}
			else
				echo "FAILURE: Autosomal coverage file ${FILE} does not exist." | tee -a "${IDENTITY}/MergeReport.txt"
			fi
			;;
		b)
			if [[ -e ${OPTARG} ]]
			then
				SRCBAM=${OPTARG}
			else
				echo "FAILURE: Incoming .BAM file ${OPTARG} does not exist." | tee -a "${IDENTITY}/MergeReport.txt"
			fi
			;;
		i)
			if [[ -e ${OPTARG} ]]
			then
				SRCBAI=${OPTARG}
			else
				echo "FAILURE: Incoming .BAI file ${OPTARG} does not exist." | tee -a "${IDENTITY}/MergeReport.txt"
			fi
			;;
		l)
			LNKBAM=${OPTARG}
			;;
		o)
			LNKBAI=${OPTARG}
			;;
		p)
			if [[ -e ${OPTARG} ]]
			then
				PLATFORM=${OPTARG}
			else
				echo "FAILURE: Platform file ${PLATFORM} does not exist." | tee -a "${IDENTITY}/MergeReport.txt"
			fi
			;;
		?)
			echo "FAILURE: ${OPTION} ${OPTARG} is not valid!" | tee -a "${IDENTITY}/MergeReport.txt"
			usage
			exit 1
			;;
	esac
done

echo "================================"  | tee -a "${IDENTITY}/MergeReport.txt"
echo "Depths of Covereage Gender comparison for ${IDENTITY}"  | tee -a "${IDENTITY}/MergeReport.txt"
echo "" | tee -a "${IDENTITY}/MergeReport.txt"

echo "" | tee -a "${IDENTITY}/MergeReport.txt"
echo "Compare X and Y chromosomal coverage as compared to an Autosomal chromosome coverage to detect chromosomal gender." | tee -a "${IDENTITY}/MergeReport.txt"
echo "If the sex chromosomes are already defined in the SampleDescriptions file, only a warning will be posted." | tee -a "${IDENTITY}/MergeReport.txt"
echo "" | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA Commandline: ${0} ${@}" | tee -a "${IDENTITY}/MergeReport.txt"
echo "" | tee -a "${IDENTITY}/MergeReport.txt"

echo "INFORMA Sample description gender: \"${GENDER}\""
echo "INFORMA Sample description chroms: \"${SEXCHR}\""
echo "INFORMA X Depth of Coverage File:  \"${XCOVER}\""
echo "INFORMA Y Depth of Coverage File:  \"${YCOVER}\""
echo "INFORMA Autosomal Coverage File:   \"${ACOVER}\""
echo "INFORMA Incoming .BAM file:        \"${SRCBAM}\""
echo "INFORMA Incoming .BAI file:        \"${SRCBAI}\""
echo "INFORMA Outgoing .BAM file:        \"${LNKBAM}\""
echo "INFORMA Outgoing .BAI file:        \"${LNKBAI}\""
echo "INFORMA Platform file:             \"${PLATFORM}\""
echo ""

if [[ -z ${XCOVER} ]] || [[ -z ${YCOVER} ]] || [[ -z ${ACOVER} ]] || [[ -z ${SRCBAM} ]] || [[ -z ${SRCBAI} ]] || [[ -z ${LNKBAM} ]] || [[ -z ${LNKBAI} ]] || [[ -z ${PLATFORM} ]]
then
	usage
	exit 1
fi

source ${PLATFORM}

Xcoverage=$(sed -n 3p ${XCOVER} | sed 's/\(^[^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)$/\3/')
Ycoverage=$(sed -n 3p ${YCOVER} | sed 's/\(^[^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)$/\3/')
Acoverage=$(sed -n 3p ${ACOVER} | sed 's/\(^[^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)\t\([^\t]\+\)$/\3/')
XARatio=$(echo "scale=3; $Xcoverage/$Acoverage" | bc | sed 's/^\./0./')
YARatio=$(echo "scale=3; $Ycoverage/$Acoverage" | bc | sed 's/^\./0./')
XCount=$(echo "scale=3; ($Xcoverage/$Acoverage)/$XRat" | bc | sed 's/^\./0./')
YCount=$(echo "scale=3; ($Ycoverage/$Acoverage)/$YRat" | bc | sed 's/^\./0./')
calculatedgender=
sexchromosomes=

echo "INFORMA Platform:   $Plat" | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA X bias:     $XRat ±$XVar" | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA Y bias:     $YRat ±$YVar" | tee -a "${IDENTITY}/MergeReport.txt"
echo ""
echo "INFORMA X coverage: $Xcoverage" | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA Y coverage: $Ycoverage" | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA A coverage: $Acoverage" | tee -a "${IDENTITY}/MergeReport.txt"
echo ""
echo "INFORMA X/A ratio:  $XARatio" | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA Y/A ratio:  $YARatio" | tee -a "${IDENTITY}/MergeReport.txt"
echo ""
echo "INFORMA X count:    $XCount" | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA Y count:    $YCount" | tee -a "${IDENTITY}/MergeReport.txt"
echo ""

XChromes=0
YChromes=0

Xmin=$(echo "scale=3;$XCount - $XVar" | bc | sed 's/^\./0./')
Xmax=$(echo "scale=3;$XCount + $XVar" | bc | sed 's/^\./0./')
Ymin=$(echo "scale=3;$YCount - $YVar" | bc | sed 's/^\./0./')
Ymax=$(echo "scale=3;$YCount + $YVar" | bc | sed 's/^\./0./')

# Count X and Y chromosomes that fall within boundries from whole numbers between 1 and 4: XXXXYYYY at most.
for i in `seq 4`
do
	XinRange=$(echo "$i < $Xmax && $i > $Xmin" | bc)
	YinRange=$(echo "$i < $Ymax && $i > $Ymin" | bc)
	
	if [ $XinRange -eq 1 ]
	then
		echo "INFORMA X:$Xmin < $XCount < $Xmax is in range of $i" | tee -a "${IDENTITY}/MergeReport.txt"
		echo ""
		XChromes=$i
	fi
	
	if [ $YinRange -eq 1 ]
	then
		echo "INFORMA Y:$Ymin < $YCount < $Ymax is in range of $i" | tee -a "${IDENTITY}/MergeReport.txt"
		echo ""
		YChromes=$i
	fi
done

# Build chromosome string.
if [ $XChromes -gt 0 ]
then
	# There are X chromosomes within the defined boundries.
	# Write a line of that many Xs.
	sexchromosomes=$(for a in `seq ${XChromes}`; do echo -n X; done)
elif [ $(echo "scale=3;$XCount > (1.0 - $XVar)" | bc) -eq 1 ]
then
	# There are no X chromosomes within the boundries.
	# Frational portions of X are greater than ONE.
	# Append these chromosomes with an E mark!
	sexchromosomes=E$(for a in `seq ${XCount}`; do echo -n X; done)
else
	# There are no X chromosomes within the boundries.
	# Frations of X found are below the lowest possible boundry.
	# Set the number of X chromoromes to ZERO.
	sexchromosomes="0"
fi

if [ $YChromes -gt 0 ]
then
	# There are Y chromosomes within the defined boundries.
	# Write a line of that many Ys
	sexchromosomes=${sexchromosomes}$(for a in `seq ${YChromes}`; do echo -n Y; done)
elif [ $(echo "scale=3;$YCount > (1.0 - $YVar)" | bc) -eq 1 ]
then
	# There are no Y chromosomes within the boundries.
	# Fraction portions of Y are greater than ONE.
	# Append these chromomes with an E mark.
	sexchromosomes=${sexchromosomes}E$(for a in `seq ${YCount}`; do echo -n Y; done)
elif [ $XChromes -eq 1 ]
then
	# There are no Y chromosomes within the boundries.
	# There is ONE X chromosome.
	# Fractional portions of Y are below lowest possible boundry.
	# Set the number of Y chromosomes to ZERO.
	sexchromosomes=${sexchromosomes}0
fi

# Decide overall gender
if [[ $XChromes -eq 0 ]]
then
	# Could not find any X chromosomes so FAIL!
	calculatedgender="Unknown"
else
	# There is at least ONE X chromosome
	if [[ $YChromes -eq 0 ]] && [[ $(echo "scale=3;$YCount < (1.0 - $YVar)" | bc) -eq 1 ]]
	then
		# There are no Y chromosomes within the boundries.
		# There are no fractional Y chromosome portions greater than the lowest possible broundry.
		calculatedgender="Female"
	else
		# There are at least 1 full Y chromosome present, even if it falls outside the boundries.
		calculatedgender="Male"
	fi
fi

echo "INFORMA Gender reported:   ${GENDER}:${SEXCHR}." | tee -a "${IDENTITY}/MergeReport.txt"
echo "INFORMA Gender calculated: $calculatedgender:$sexchromosomes." | tee -a "${IDENTITY}/MergeReport.txt"
echo ""


if [[ ${SEXCHR} != "" ]] && [[ ${SEXCHR} != $sexchromosomes ]]
then
	# Sex chromosomes where specified
	# Specified sex chromosomes do not match calculated sex chromosomes.
	# Don't FAIL but WARN.
	echo "WARNING Calculated chromosomal gender $sexchromosomes does not match recorded chromosomal gender ${SEXCHR}" | tee -a "${IDENTITY}/MergeReport.txt"
	echo "WARNING The sample will be processed according to recorded chromosomal gender ${SEXCHR}" | tee -a "${IDENTITY}/MergeReport.txt"
elif [ $calculatedgender == "Unknown" ]
then
	# Unable to determine gender. No X or Ys detected.
	# FAILS.
	echo "FAILURE Possibly wrong capture platform specified, poor capture, contamination of sample or mosaic aneuploidy" | tee -a "${IDENTITY}/MergeReport.txt"
	exit 1
elif [[ $GENDER == "Unknown" ]] && [[ $GENDER != $calculatedgender ]]
then
	# Incoming gender was UNKNOWN and calculated gender is KNOWN.
	# FAILS, but in a good way!
	echo "INFORMA Gender now identified!" | tee -a "${IDENTITY}/MergeReport.txt"
	echo "INFORMA Please insert this data into the Sample Description file under either Gender or Sex Chromosomes, or both." | tee -a "${IDENTITY}/MergeReport.txt"
	echo "INFORMA You can simply re-launch the pipeline to pick up where you left off." | tee -a "${IDENTITY}/MergeReport.txt"
	exit 1
elif [[ ${SEXCHR} == "UNKNOWN" ]] && [[ $sexchromosomes != "" ]] && [[ $sexchromosomes != *"E"* ]]
then
	# Sex chromosomes where unknown (majority)
	# Calculated sex chromosomes are detected
	# Calculated sex chromosomes are all within boundries. No E in sequence.
	echo "INFORMA Chromosomal gender calculated!" | tee -a "${IDENTITY}/MergeReport.txt"
	echo "INFORMA You can specify this in the Sample description file." | tee -a "${IDENTITY}/MergeReport.txt"
elif [[ $sexchromosomes != "XX" ]] && [[ $sexchromosomes != "XY" ]] # && [[ $sexchromosomes != "" ]]
then
	# Calculated Sex chromosomes something other than XX or XY
	# FAIL.
	echo "FAILURE Possible sex chromosome aneuploidy." | tee -a "${IDENTITY}/MergeReport.txt"
	echo "FAILURE If you wish to proceed, please enter this or your preferred gender into the Sample Descriptions file's Sex Chromosomes column, then re-run this pipeline." | tee -a "${IDENTITY}/MergeReport.txt"
	exit 1
elif [[ ${SEXCHR} == $sexchromosomes ]]
then
	# Sex chromosomes where specified and match calculated sex chromosomes.
	echo "SUCCESS Calculated chromosomal gender $sexchromosomes matches recorded chromosomal gender ${SEXCHR}" | tee -a "${IDENTITY}/MergeReport.txt"
fi

if [[ $GENDER != $calculatedgender ]]
then
	# Specified gender doesn't match calculated gender.
	# FAIL.
	echo "FAILURE Calculated gender $calculatedgender conflicts with reported gender $GENDER" | tee -a "${IDENTITY}/MergeReport.txt"
	echo "FAILURE Please verify your sample authenticity or correct the Sample Desciption file." | tee -a "${IDENTITY}/MergeReport.txt"
	exit 1
else
	echo "SUCCESS Calculated gender matches reported gender so processing normally!" | tee -a "${IDENTITY}/MergeReport.txt"
fi

#exit 1

${HARDLINKER} ${SRCBAM} ${LNKBAM} && ${HARDLINKER} ${SRCBAI} ${LNKBAI} && exit 0

exit 1

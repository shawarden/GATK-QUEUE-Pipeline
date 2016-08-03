#!/bin/bash

SOURCE=${1}
LINKED=${2}

LINKED_DIR=$(dirname ${LINKED})
if mkdir -p ${LINKED_DIR}
then	# Created link destination, just in case it's different.
	if [ -f ${SOURCE} ]
	then	# source file exists!
		if [ -f ${LINKED} ]
		then	# Link target file exists!
			if [ ${LINKED} -ef ${SOURCE} ]
			then	# Linked file is already mapped to the same inode as the Source.
				echo -e "SUCCESS: ${LINKED} is already hard linked to ${SOURCE}\n"
				exit 0
			else	# Linked file is not mapped to the same inode as the Source.
				echo -e "WARNING: ${LINKED} already exists but is not linked to ${SOURCE}. Deleting...\n"
				if rm --interactive=never -r ${LINKED}
				then	# successfully deleted original link file.
					if ${0} ${@}	# recursion ftw!
					then
						exit 0
					fi
				else	# could not delete original link file.
					echo -e "FAILURE: Unable to delete ${LINKED} file. Check permissions.\n"
				fi
			fi
		else	# Link target file doesn't exist.
			if ln ${SOURCE} ${LINKED}
			then	# created hard link!
				echo -e "SUCCESS: ${SOURCE} hard linked to ${LINKED}.\n"
				exit 0
			else	# Unable to create hard link.
				echo -e "FAILURE: Unable to hard linked ${SOURCE} to ${LINKED}. Check permissions. Linked cannot be made across file-systems.\n"
			fi
		fi
	else	# source file doesn't exist?
		echo -e "FAILURE: Source file ${SOURCE} doesn't exist!\n"
	fi
else	# unable to create link destination folder.
	echo -e "FAILURE: Unable to create ${LINKED} folder: ${LINKED_DIR}. Please check folder permissions.\n"
fi

exit 1	# if we got here, something's wrong.
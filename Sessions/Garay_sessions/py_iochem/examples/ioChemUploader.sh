#!/bin/bash

# User parameters: take from bash script config file
readarray -t iochem_info < "bash_upload_config.dat"
remote_folder=${iochem_info[0]}
remote_desc=${iochem_info[1]}
shell_route=${iochem_info[2]}
iochem_route=${iochem_info[3]}

# Read array with calc information: name, input and output file
readarray -t upinfo < "iochem_upload_summ.dat"

# ioChem connection
source ${shell_route}/start-rep-shell
exe-rep-command cdpro ${iochem_route}
exe-rep-command pwdpro
exe-rep-command cpro -n ${remote_folder} -d ${remote_desc}
exe-rep-command cdpro ${remote_folder}
exe-rep-command pwdpro

initfold=$(pwd)
echo "" > calcid_tracking.dat
mkdir -p upload_information
# Extract info
for line in "${upinfo[@]}" ; do
	echo "=> $line"
	IFS=";" read -a arr <<< $line
	name=${arr[0]}
	ifile=${arr[1]}
	ofile=${arr[2]}
	desc="${name}"
	if [[ -f ${ifile} ]] && [[ -f ${ofile} ]] ; then
		echo "loadcalc -n ${name} -d ${desc} -i $ifile -o $ofile"
		exe-rep-command loadcalc -n ${name} -d ${desc} -i $ifile -o $ofile > upload_information/${name}_upload_info.dat
		# get the ID
		calcid=$(tail -n1 upload_information/${name}_upload_info.dat)
		echo "${name} ; ${calcid}" >> calcid_tracking.dat
	else
		echo "Did not find $ifile for $name"
	fi
done
# cleanup
tail -n+2 calcid_tracking.dat > calcid_tracking.temp
mv calcid_tracking.temp calcid_tracking.dat
# go back to parent folder
# remote (not required, but keep for safety)
exe-rep-command cdpro "${iochem_route}/${remote_folder}"
# local
cd ${initfold}



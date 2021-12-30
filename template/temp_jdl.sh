#!/bin/sh

#define parameters which are passed in.
jdlmain_p=$1
jdlstor_p=$2
jdlparticle=$3
jdlenergy=$4
jdlevent=$5
jdltag=$6
jdlcfg=$7
jdlmac=$8
#define the template.
#on_exit_hold =   (ExitCode != 0)
#on_exit_remove =   (ExitCode == 0)
#job_machine_attrs = Machine
#job_machine_attrs_history_length = 3

#Output = ${jdlpath}/run_${jdlenergy}GeV_${jdlparticle}_N${jdlevent}_${jdltag}.out
#Error =  ${jdlpath}/run_${jdlenergy}GeV_${jdlparticle}_N${jdlevent}_${jdltag}.err
#Log =    ${jdlpath}/run_${jdlenergy}GeV_${jdlparticle}_N${jdlevent}_${jdltag}.log

#Output = ${jdlpath}/condor_output/sce_\$(cluster)_\$(process).stdout
#Error =  ${jdlpath}/condor_output/sce_\$(cluster)_\$(process).stderr
#Log =    ${jdlpath}/condor_output/sce_\$(cluster)_\$(process).condor
#Requirements = TARGET.Machine =?= "r510-0-1.privnet"
#Requirements = TARGET.Machine =?= "r510-0-1.privnet" ||  TARGET.Machine =?= "r510-0-9.privnet" || TARGET.Machine =?= "r720-0-2.privnet" || TARGET.Machine =?= "r720-datanfs.privnet"
#Requirements = TARGET.Machine =?= "r510-0-1.privnet" ||  TARGET.Machine =?= "r510-0-9.privnet" || TARGET.Machine =?= "r720-0-2.privnet" || TARGET.Machine =?= "r720-datanfs.privnet"
cat  << EOF
universe = vanilla
getenv=True
Executable = ${jdlmain_p}/CEPC_CaloTiming 
should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
#Requirements = TARGET.Machine =?= "r510-0-1.privnet" ||  TARGET.Machine =?= "r510-0-9.privnet" || TARGET.Machine =?= "r720-0-2.privnet" || TARGET.Machine =?= "r720-datanfs.privnet"
Output = ${jdlstor_p}/run_${jdlenergy}GeV_${jdlparticle}_N${jdlevent}_${jdltag}.out
Error =  ${jdlstor_p}/run_${jdlenergy}GeV_${jdlparticle}_N${jdlevent}_${jdltag}.err
Log =    ${jdlstor_p}/run_${jdlenergy}GeV_${jdlparticle}_N${jdlevent}_${jdltag}.log
Arguments = -c ${jdlcfg} -m ${jdlmac} -o ${jdlstor_p}/${jdlenergy}GeV_N${jdlevent}_${jdlParticle}_${jdltag}
Queue 1
EOF

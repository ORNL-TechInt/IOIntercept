*** Settings ***
Library           spectralRobotLibrary.py	${WORKDIR}
Library           OperatingSystem
Library           Collections
Library		  Process
Suite Setup       Start Spectral Unit Tests

*** Variables ***
${WORKDIR}	/tmp/spectral
${EXDIR}	${WORKDIR}/examples
${GTCDIR}	${WORKDIR}/examples/gtc
${FILECNT}	5
${MPIPATH}	/lib64/openmpi
${MPIBIN}	/lib64/openmpi/bin

*** Test Cases ***
Build Spectral
    [Documentation]	Build Spectral 
    [Tags]	Setup 
    Build Spectral
    Build Outcome

Test Simple Static
    [Documentation]	Test Static
    [Tags]	Simple
    Set Environment Variable	PFS_DIR	/tmp/gpfs
    Set Environment Variable	PERSIST_DIR	/mnt/bblv_57
    Clean Directory	/tmp/gpfs
    Create Directory	/tmp/gpfs
    Run Simple	Static	${EXDIR}	${FILECNT}
    Status Should Be	${FILECNT}

Test Simple Dynamic
    [Documentation]	Test Dynmaic
    [Tags]	Simple
    Set Environment Variable	PFS_DIR	/tmp/gpfs
    Set Environment Variable	PERSIST_DIR	/mnt/bblv_57
    Clean Directory	/tmp/gpfs
    Create Directory	/tmp/gpfs
    Run Simple	Dynamic	${EXDIR}	${FILECNT}
    Status Should Be	${FILECNT}

Test MPI Static
    [Documentation]	Simple Static MPI Test
    [Tags]	MPI
    Set Environment Variable	PFS_DIR	/tmp/gpfs
    Set Environment Variable	PERSIST_DIR	/mnt/bblv_57
    Clean Directory	/tmp/gpfs
    Create Directory	/tmp/gpfs
    Run MPI	${EXDIR}	${MPIBIN}
    Status Should Be	4

Build GTC
    [Documentation]	Build GTC
    [Tags]	GTC
    Build GTC	${GTCDIR}
    Build Outcome

*** Keywords *** 
Build Spectral
    [Arguments]
    [Timeout]	10 seconds
    Build System	${WORKDIR}

Build Outcome
    [Arguments]
    [Timeout]	10 seconds
    Status Should Not Contain	error
    Status Should Not Contain	failure

Create Directory
    [Arguments]	${dest_dir}
    [Timeout]	10 seconds
    Run Process	mkdir	-p	${dest_dir}

Clean Directory
    [Arguments]	${dest_dir}
    [Timeout]	10 seconds
    Run Process	rm	-rf	${dest_dir}

Start Spectral Unit Tests
    log	Starting Spectral Unit Testing

End Spectral Unit Tests
    log	Terminating Spectral Unit Testing

*** Settings ***
Suite Teardown    End Spectral Unit Tests

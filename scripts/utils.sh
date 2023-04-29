#!/bin/bash

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OUTPUT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Use colored stdout if the terminal supports it
## and as long as a stdout are not redirected
function term_colors {
    local action
    local what
    
    action=$1
    what=$2
    [[ -z "${what}" ]] && what=1
    
    if [[ "${action}" == "enable" && -t "${what}" ]]; then
        ## ANSI foreground colors
        black=$(tput setaf 0)
        red=$(tput setaf 1)
        green=$(tput setaf 2)
        yellow=$(tput setaf 3)
        blue=$(tput setaf 4)
        magenta=$(tput setaf 5)
        cyan=$(tput setaf 6)
        white=$(tput setaf 7)

        ## Text modes
        bold=$(tput bold)
        dim=$(tput dim)
        reset=$(tput sgr0)
    else
        export black=
        export red=
        export green=
        export yellow=
        export blue=
        export magenta=
        export cyan=
        export white=

        export bold=
        export dim=

        export reset=
    fi
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TESTING
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function test_context {
    local magenta
    local reset
    if [[ -t 1 ]]; then
        magenta=$(tput setaf 5)
        reset=$(tput sgr0)
    fi
    echo -e "${magenta}*** $*${reset}"
}

function test_context_begin {
    test_context "$* ..."
}

function test_context_end {
    test_context "$* ... done"
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CONDITIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function error {
    local bold
    local gray
    local red
    local reset
    
    ON_ERROR=${ON_ERROR:-on_error}
    TRACEBACK_ON_ERROR=${TRACEBACK_ON_ERROR:-true}
    EXIT_ON_ERROR=${EXIT_ON_ERROR:-true}
    EXIT_VALUE=${EXIT_VALUE:-1}

    ## Parse arguments
    while [ -n "$1" ]; do
        case "$1" in
            --dryrun) EXIT_ON_ERROR=false; shift;;
            --value=*) EXIT_VALUE="${1/--value=/}"; shift;;
            *) break;;
        esac
    done

    if [[ -t 1 ]]; then
        red=$(tput setaf 1)
        gray=$(tput setaf 8)
        bold=$(tput bold)
        reset=$(tput sgr0)
    fi

    echo -e "${red}${bold}ERROR:${reset} ${bold}$*${reset}"

    if ${TRACEBACK_ON_ERROR}; then
       echo -e "${gray}Traceback:"
       for ((ii = 1; ii < "${#BASH_LINENO[@]}"; ii++ )); do
           printf "%d: %s() on line #%s in %s\\n" "$ii" "${FUNCNAME[$ii]}" "${BASH_LINENO[$((ii-1))]}" "${BASH_SOURCE[$ii]}"
       done
    fi

    if [[ -n "${ON_ERROR}" ]]; then
        if [[ $(type -t "${ON_ERROR}") == "function" ]]; then
            ${ON_ERROR}
        fi
    fi

    ## Exit?
    if ${EXIT_ON_ERROR}; then
        echo -e "Exiting (exit ${EXIT_VALUE})${reset}";
        exit "${EXIT_VALUE}"
    fi

    printf "%s" "${reset}"
}

function warn {  
    local bold
    local yellow
    local reset
    
    TRACEBACK_ON_WARN=${TRACEBACK_ON_WARN:-false}
    
    if [[ -t 1 ]]; then
        yellow=$(tput setaf 3)
        bold=$(tput bold)
        reset=$(tput sgr0)
    fi
    
    echo -e "${yellow}${bold}WARNING${reset}: $*"
    
    if ${TRACEBACK_ON_WARN}; then
       echo -e "${gray}Traceback:"
       for ((ii = 1; ii < "${#BASH_LINENO[@]}"; ii++ )); do
           printf "%d: %s() on line #%s in %s\\n" "$ii" "${FUNCNAME[$ii]}" "${BASH_LINENO[$((ii-1))]}" "${BASH_SOURCE[$ii]}"
       done
    fi
    
    printf "%s" "${reset}"
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ASSERTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Usage: assert_file_exists /path/to/file
function assert_file_exists {
    local file
    local name
    local prefix
    
    [[ $# -ne 1 ]] && error "${FUNCNAME[0]}() requires a single argument: $#"
    file=$1
    prefix=""
    if [[ $file == --by-env=* ]]; then
	name=${file##--by-env=}
	file=${!name}
	prefix="$name="
    fi
    [[ -n "$file" ]] || error "File name must be non-empty: $prefix'$file'"
    [[ -f "$file" ]] || error "No such file: $prefix'$file' (working directory '${PWD}')"
}

function assert_link_exists {
    local file
    local name
    local prefix
    
    [[ $# -ne 1 ]] && error "${FUNCNAME[0]}() requires a single argument: $#"
    [[ -n "$1" ]] || error "File name must be non-empty: '$1'"
    file=$1
    prefix=""
    if [[ $file == --by-env=* ]]; then
	name=${file##--by-env=}
	file=${!name}
	prefix="$name="
    fi
    [[ -L "$file" ]] || error "File is not a link: $prefix'$file' (working directory '${PWD}')"
    [[ -e "$file" ]] || error "[File] link is broken: $prefix'$file' (working directory '${PWD}')"
}

## Usage: assert_file_executable /path/to/file
function assert_file_executable {
    local file
    local name
    local prefix
    
    [[ $# -ne 1 ]] && error "${FUNCNAME[0]}() requires a single argument: $#"
    assert_file_exists "$1"
    file=$1
    prefix=""
    if [[ $file == --by-env=* ]]; then
	name=${file##--by-env=}
	file=${!name}
	prefix="$name="
    fi
    [[ -x "$file" ]] || error "File exists but is not executable: $prefix'$file' (working directory '${PWD}')"
}

## Usage: assert_directory_exists /path/to/folder
function assert_directory_exists {
    local dir
    local name
    local prefix
    
    [[ $# -ne 1 ]] && error "${FUNCNAME[0]}() requires a single argument: $#"
    dir=$1
    prefix=""
    if [[ $dir == --by-env=* ]]; then
	name=${dir##--by-env=}
	dir=${!name}
	prefix="$name="
    fi
    [[ -n "$dir" ]] || error "Directory name must be non-empty: $prefix='$dir'"
    [[ -d "$dir" ]] || error "No such directory: $prefix='$dir' (working directory '${PWD}')"
}

## Usage: assert_patient_name "${PATIENT}"
function assert_patient_name {
    [[ $# -ne 1 ]] && error "${FUNCNAME[0]}() requires a single argument: $#"
    [[ -n "$1" ]] || error "Patient name must be non-empty: '$1'"
    [[ "$1" == *[_]* ]] && error "Patient name must not contain underscores: $1"
}


## Usage: assert_pwd
function assert_pwd {
    ## Need to disable this assertion to avoid triggering a false SC2119 ShellCheck warning
    ## [[ $# -ne 0 ]] && error "${FUNCNAME[0]}() must not be called with arguments: $#"
    
    ## Don't allow running the pipeline from within LG3_HOME
    equal_dirs "${PWD}" "${LG3_HOME}" && error "The LG3 Pipeline must not be run from the folder where it is installed (LG3_HOME): ${PWD}"
}


## Usage: assert_python "" or assert_python "<python-binary>"
function assert_python {
    local bin
    local name
    local prefix
    local version
    local version_x_y
    
    ## Arguments are optional
    bin=$1

    prefix=""
    if [[ $bin == --by-env=* ]]; then
	name=${bin##--by-env=}
	bin=${!name}
	prefix="$name="
    fi
    
    
    if [[ -n "$bin" ]]; then
        assert_file_executable "$1"
    else
        bin=$(command -v python) || error "Python executable not found on PATH: ${PATH}"
    fi
    
    ## Assert correct version
    version=$(2>&1 "$bin" --version | sed -E 's/.*(P|p)ython *//g')
    version_x_y=$(echo "$version" | sed -E 's/[.][0-9]+$//g')
    [[ "$version_x_y" == "2.6" ]] || [[ "$version_x_y" == "2.7" ]] || error "Requires Python 2.6 or 2.7: $version ($bin)"
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NAVIGATION
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function change_dir {
    local opwd
    
    opwd=${PWD}
    assert_directory_exists "$1"
    cd "$1" || error "Failed to set working directory to $1"
    echo "New working directory: '$1' (was '${opwd}')"
}

function make_dir {
    mkdir -p "$1" || error "Failed to create new working directory $1"
}

function make_change_dir {
    make_dir "$1"
    change_dir "$1"
}

function equal_dirs {
    local a
    local b
    a=$(readlink -f "$1")
    b=$(readlink -f "$2")
    [[ "${a}" == "${b}" ]]
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LG3 specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function source_lg3_conf {
    ## Already sourced?
    if [[ -n ${LG3_CONF_SOURCED} ]]; then
	echo "All lg3.conf files have already been sourced and validates. Skipping. [LG3_CONF_USER='${LG3_CONF_USER}']"
	return 0;
    fi
	
    ## The default settings
    assert_file_exists "${LG3_HOME}/lg3.conf"
    # shellcheck disable=1090
    source "${LG3_HOME}/lg3.conf"
    echo "Sourced: ${LG3_HOME}/lg3.conf"

    ## lg3.conf per LG3_CONF_USER?
    if [[ -z ${LG3_CONF_USER} ]]; then
        if [ -f "lg3.conf" ] && ! equal_dirs "." "${LG3_HOME}"; then
            echo "Located: ./lg3.conf [PWD=${PWD}]"
            LG3_CONF_USER=${PWD}/lg3.conf
	fi
    else
        echo "Using: LG3_CONF_USER='${LG3_CONF_USER}'"
    fi
    
    if [[ -n ${LG3_CONF_USER} ]]; then
	assert_file_exists "${LG3_CONF_USER}"
	
        # shellcheck disable=1090
        source "${LG3_CONF_USER}"
        echo "Sourced: ${LG3_CONF_USER} ($(stat --printf='%s' "${LG3_CONF_USER}") bytes)"
    fi

    ## Assert that all LG3 software dependencies are valid
    lg3_assert_software
    
    LG3_CONF_SOURCED=true
    export LG3_CONF_USER
}



function lg3_software_envvars {
    echo "JAVA PYTHON RSCRIPT ANNOVAR_HOME BEDTOOLS BWA CUTADAPT GATK MUTECT PICARD_HOME SAMTOOLS PINDEL PINDEL2VCF TG"
}

function lg3_list_software {
    local name
    local value
    
    echo "LG3 software dependencies used:"
    for name in $(lg3_software_envvars); do
        value=${!name}
	if [[ -z "${value}" ]]; then
            value="<empty>"
	else
            value="'${value}'"
	fi
        echo "- ${name}=${value}"
    done
}

function lg3_assert_software {
    assert_file_executable  --by-env=JAVA
    assert_file_executable  --by-env=PYTHON
    assert_file_executable  --by-env=PYTHON
    assert_python           --by-env=PYTHON
    assert_file_executable  --by-env=RSCRIPT
    assert_directory_exists --by-env=ANNOVAR_HOME
    assert_file_executable  --by-env=BEDTOOLS
    assert_file_executable  --by-env=BWA
    assert_file_executable  --by-env=CUTADAPT
    assert_file_exists      --by-env=GATK
    assert_file_exists      --by-env=MUTECT
    assert_directory_exists --by-env=PICARD_HOME
    assert_file_executable  --by-env=SAMTOOLS
    assert_file_executable  --by-env=PINDEL
    assert_file_executable  --by-env=PINDEL2VCF
    assert_file_executable  --by-env=TG
}

function lg3_qsub_envvar_append_software {
    local name
    local value
    local res

    lg3_assert_software
    
    res="${QSUB_ENVVARS}"
    for name in $(lg3_software_envvars); do
        value=${!name}
	if [[ -n "${value}" ]]; then
            res="${res},${name}=${value}"
	fi
    done
    
    echo "${res}"
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Job scheduler specifics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Does 'qsub' on the current system support CLI option -d <path>?
qsub_can_set_pwd=false
if [[ ${CLUSTER} == "tipcc" ]]; then
    # shellcheck disable=2034
    qsub_can_set_pwd=true
fi


if [[ -n "${PBS_JOBID}" ]]; then
    LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:-/scratch/${USER:?}/${PBS_JOBID}}
elif [[ -n "${SLURM_JOB_ID}" ]]; then
    LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:-/scratch/${USER:?}/${SLURM_JOB_ID}}
elif [[ -n "${TMPDIR}" ]]; then
    LG3_SCRATCH_ROOT=${LG3_SCRATCH_ROOT:-${TMPDIR}}
else
    error "Cannot set LG3_SCRATCH_ROOT"
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## FIXME: Can we remove this?
unset PYTHONPATH  ### ADHOC: In case it is set by user. /HB 2018-09-07

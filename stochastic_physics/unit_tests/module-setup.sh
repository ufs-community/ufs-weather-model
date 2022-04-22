# Create a test function for sh vs. bash detection.  The name is
# randomly generated to reduce the chances of name collision.
__ms_function_name="setup__test_function__$$"
eval "$__ms_function_name() { /bin/true ; }"

# Determine which shell we are using
__ms_ksh_test=$( eval '__text="text" ; if [[ $__text =~ ^(t).* ]] ; then printf "%s" ${.sh.match[1]} ; fi' 2> /dev/null | cat )
__ms_bash_test=$( eval 'if ( set | grep '$__ms_function_name' | grep -v name > /dev/null 2>&1 ) ; then echo t ; fi ' 2> /dev/null | cat )

if [[ ! -z "$__ms_ksh_test" ]] ; then
    __ms_shell=ksh
elif [[ ! -z "$__ms_bash_test" ]] ; then
    __ms_shell=bash
else
    # Not bash or ksh, so assume sh.
    __ms_shell=sh
fi

if [[ -d /lfs3 ]] ; then
    # We are on NOAA Jet
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    module purge
elif [[ -d /scratch1 && ! -d /scratch ]] ; then
    # We are on NOAA Hera
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    module purge
elif [[ -d /work/noaa ]] ; then
    # We are on Orion
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /apps/lmod/init/$__ms_shell
    fi
    module purge
elif [[ -d /data ]] ; then
    # We are on SSEC Wisconsin S4
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /usr/share/lmod/lmod/init/$__ms_shell
    fi
    module purge
elif [[ -d /gpfs/hps && -e /etc/SuSE-release ]] ; then
    # We are on NOAA Luna or Surge
    if ( ! eval module help > /dev/null 2>&1 ) ; then
       source /opt/modules/default/init/$__ms_shell
    fi
    module purge
    module purge
    # Workaround until module issues are fixed:
    unset _LMFILES_
    unset LOADEDMODULES
    module use /opt/modulefiles
    module use /opt/cray/ari/modulefiles
    module use /opt/cray/craype/default/alt-modulefiles
    module use /opt/cray/alt-modulefiles
    module use /gpfs/hps/nco/ops/nwprod/modulefiles
    module use /gpfs/hps/nco/ops/nwprod/lib/modulefiles
    module use /usrx/local/prod/modulefiles
elif [[ -d /dcom && -d /hwrf ]] ; then
    # We are on NOAA Tide or Gyre
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /usrx/local/Modules/default/init/$__ms_shell
    fi
    module purge
elif [[ -L /usrx && "$( readlink /usrx 2> /dev/null )" =~ dell ]] ; then
    # We are on NOAA Mars or Venus
    if ( ! eval module help > /dev/null 2>&1 ) ; then
	source /usrx/local/prod/lmod/lmod/init/$__ms_shell
    fi
    module purge
elif [[ -d /glade ]] ; then
    # We are on NCAR Cheyenne
    if ( ! eval module help > /dev/null 2>&1 ) ; then
      . /glade/u/apps/ch/modulefiles/default/localinit/localinit.sh
    fi
    module purge
elif [[ -d /work/stampede ]] ; then
    # We are on TACC Stampede
    if ( ! eval module help > /dev/null 2>&1 ) ; then
      source /opt/apps/lmod/lmod/init/bash
    fi
    module purge
elif [[ -d /lustre && -d /ncrc ]] ; then
    # We are on GAEA. 
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        # We cannot simply load the module command.  The GAEA
        # /etc/profile modifies a number of module-related variables
        # before loading the module command.  Without those variables,
        # the module command fails.  Hence we actually have to source
        # /etc/profile here.
        source /etc/profile
        __ms_source_etc_profile=yes
    else
        __ms_source_etc_profile=no
    fi
    module purge
# clean up after purge 
    unset _LMFILES_
    unset _LMFILES_000
    unset _LMFILES_001
    unset LOADEDMODULES
    module load modules
    if [[ -d /opt/cray/ari/modulefiles ]] ; then
        module use -a /opt/cray/ari/modulefiles
    fi
    if [[ -d /opt/cray/pe/ari/modulefiles ]] ; then
        module use -a /opt/cray/pe/ari/modulefiles
    fi
    if [[ -d /opt/cray/pe/craype/default/modulefiles ]] ; then
        module use -a /opt/cray/pe/craype/default/modulefiles
    fi
    if [[ -s /etc/opt/cray/pe/admin-pe/site-config ]] ; then
        source /etc/opt/cray/pe/admin-pe/site-config
    fi
    export NCEPLIBS=/lustre/f1/pdata/ncep_shared/NCEPLIBS/lib
    if [[ -d "$NCEPLIBS" ]] ; then
        module use $NCEPLIBS/modulefiles
    fi
    if [[ "$__ms_source_etc_profile" == yes ]] ; then
      source /etc/profile
      unset __ms_source_etc_profile
    fi
elif [[ -d /Applications ]] ; then
    # We are on a MacOSX system, nothing to do
    echo "Platform: MacOSX"
elif [[ -e /etc/redhat-release ]] ; then
    if grep -iq centos "/etc/redhat-release" ; then
        # We are on CentOS Linux, nothing to do
        echo "Platform: CentOS Linux"
    else
        echo WARNING: UNKNOWN PLATFORM 1>&2
    fi
elif [[ -e /etc/issue ]] ; then
    if grep -iq ubuntu "/etc/issue" ; then
        # We are on Ubuntu Linux, nothing to do
        echo "Platform: Ubuntu Linux"
    else
        echo WARNING: UNKNOWN PLATFORM 1>&2
    fi
else
    echo WARNING: UNKNOWN PLATFORM 1>&2
fi

unset __ms_shell
unset __ms_ksh_test
unset __ms_bash_test
unset $__ms_function_name
unset __ms_function_name

.. _UsingGEOSOnSherlock:

Using GEOS on Stanford *Sherlock* cluster
===========================================

*Sherlock*  is the new Stanford mixed CPU/GPU `cluster<https://www.sherlock.stanford.edu/docs/tech/#in-a-nutshell>`_ that is usable on request
for students. 

There is two main user profile on *Sherlock* that can be interested in different
deployement of *GEOS*:

    - _Users_, who are only interested in running the current version of develop
    on multi-core/multi-gpu arch to get the most of the HPC.
    - _Developers_ and code wonderers, who wants to explore code, run test branches
    and even develop on premise. 

Both profile will be interested in finding pre-compiled and deployed version of
the most recent thirdPartyLibs (to the exception of development in these thirdPartyLibs).

Auto-deployement of pre-compiled binaries
---------------------------------------------

To that extend and also to fit specific needs of _Users_ profile, an automated deployement
is done thanks the following script:

.. code-block::
    :caption: Auto-deployement of GEOS on Sherlock

    #!/bin/sh

    ## 1. get commit hash and dl bundle
    # for now as it is travis commit only
    # after merge you can pull Sherlock bundle directly
    export commit_hash=$(curl -fsSL -H "Accept: application/vnd.github+json" https://api.github.com/repos/GEOS-DEV/GEOS/branches/develop | jq -r ".commit.sha")
    # or force commit
    export commit_hash=611031a


    ## make sure the commit is older that a day
    export commit_day=$(curl -fsSL -H "Accept: application/vnd.github+json" https://api.github.com/repos/GEOS-DEV/GEOS/commits/${commit_hash} | jq -r '.commit.committer.date')
    #if commit younger than one day it did not yet pass the nightly test hence we have to find the parents until it is true
    while [ $(date -d $(date -I) +%s) -le $(date -d ${commit_day:0:-1} +%s) ]
    do

        echo "commit ... ${commit_hash:0:7} ... ${commit_day}"
        export commit_hash=$(curl -fsSL -H "Accept: application/vnd.github+json" https://api.github.com/repos/GEOS-DEV/GEOS/commits/${commit_hash} | jq -r '.parents[0].sha')
        export commit_day=$(curl -fsSL -H "Accept: application/vnd.github+json" https://api.github.com/repos/GEOS-DEV/GEOS/commits/${commit_hash} | jq -r '.commit.committer.date')

    done


    echo "commit ... ${commit_hash:0:7} ... ${commit_day}"

    # loop on devices
    for DEVICE in CPU 
        do

        cd ${DEVICE}

        #getting the bundle and untar inplace
        RET=$(curl -sIL  -w"%{http_code}" https://storage.googleapis.com/geosx/Sherlock-${DEVICE}/GEOSX-and-TPL-${commit_hash:0:7}.tar.gz -o /dev/null) 
        echo "return curl .. ${RET}"
            if [[ ${RET} -eq 200 ]] 
            then
                curl -L https://storage.googleapis.com/geosx/Sherlock-${DEVICE}/GEOSX-and-TPL-${commit_hash:0:7}.tar.gz | tar --strip-components=1 --keep-old-files -xzf - 
            else
                echo "Failed. Step : download\n"
                exit 1
            fi
        ## 2. list linkage after module load and play a small test to be sure it is working
        module purge -f
        module restore geosx-cpu

        export GEOSX_EXE=${PWD}/GEOSX-${commit_hash:0:7}/bin/geosx

        ldd ${GEOSX_EXE} >> log.GEOSX-${commit_hash:0:7}.linkage
            if grep "found" log.GEOSX-${commit_hash:0:7}.linkage; then
                echo "--  FAILED LINKAGE" >>  log.GEOSX-${commit_hash:0:7}.linkage
            fi


        ERR=$(${GEOSX_EXE} --help)

            if [[ $ERR -gt 0 ]]
            then
                "GEOSX failed"
            fi


        done

    exit 0

This will pull from a google-cloud storage the correct version that have been uploaded
there as a part of the CI workflow that moves the resulting artifacts to that cloud location.

The path to which to find these deployed versions is `/oak/stanford/schools/ees/COLLABORATIONS/geosx/CPU`.
The folders name norm is `GEOSX_TPL-${TAG}-${commit_hash}` for TPL and `GEOS-${geos_commit_hash}` for GEOS. 
Indeed, TPL versions changing on a least frequent basis as GEOS, the TPL version is tested at unpack and not
inflated is not necessary.


.. note::
    `/oak/stanford/schools/ees/COLLABORATIONS/geosx/GPU/` will be recieving the GPU
    deployement as soon as it is available.


Using singularity
-------------------

An other possibility for _Users_ is to use containerized version of the code
that are made available on Sherlock through a simple succession of few commands as
described in the :ref:`UsingSingularity`.

Local version
--------------

This last option is more directed to _Developers_. Following the cloning steps from
:ref:`QuickStart`. 

A version of such a deployement can be found at `/oak/stanford/schools/ees/COLLABORATIONS/geosx/compiled`
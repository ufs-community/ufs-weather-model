Introduction
============

The three branches dtc/hwrf-physics in ccpp-physics, fv3atm and ufs-weather-model
were created from the hafs-community repository on April 3, 2020:

git clone -b support/HAFS --recursive https://github.com/hafs-community/ufs-weather-model

The original hashes were:

ccpp-physics      31a294a610a538ade4670d9e32856c4a2b7e59ad
fv3atm            23c49d0df3dee47bd7d13f0cbefc5de23c7b2663
ufs-weather-model c2643caf49ff350b66f35fd72afbcdaa65f609bd


PRs pulled into the dtc/hwrf-physics branches
=============================================

The following PRs were pulled into these branches on the given dates
and with the given hashes (further down is later). Additional commits
were made by @climbfuji as required for changing the target branch
(from dtc/develop to dtc/hwrf-physics), bug fixes, cleanup, and as
a result of the code review processes.

ccpp-physics
------------

04/06/2020 MONINEDMF https://github.com/NCAR/ccpp-physics/pull/395 -> https://github.com/NCAR/ccpp-physics/pull/428
    git cherry-pick 942889ab0acace519254712cb4e14b6aaf3e0415
    git cherry-pick 500d53a21027f362a9c12c10767f0b4f8cf3361c
    git cherry-pick b0f04b210bd588673182023ea36b56ee94642c3e
    git cherry-pick e895c62cf49b908927678d72c1ad21ad828e3587

04/08/2020 RRTMG https://github.com/NCAR/ccpp-physics/pull/412 -> https://github.com/NCAR/ccpp-physics/pull/430
    git cherry-pick 9e9222a470dd8a644ba4c3d010766c62dc71de59
    git cherry-pick 4d9e68fc27bb3974f1a6c5733bf1b8127165581d
    git cherry-pick 6f9fec9d183531c81475b3496504dc3225011472
    git cherry-pick 57873f2045e1cd2d09c71481c0593f3d7eee24ef
    git cherry-pick 7e492ca69cab20d2457af82c5d535d8b64e2ab1e
    git cherry-pick 298d1aed7b2f29ccf624a71be7d9f804eaf9e8ab
    git cherry-pick 5404462a72fe10477595c25baab0ae28fe667f0f

fv3atm
------

04/06/2020 MONINEDMF https://github.com/NCAR/fv3atm/pull/18 -> https://github.com/NCAR/fv3atm/pull/37
    git cherry-pick 80ce8591b6af615c9f85a2b3becf8aadc836d734

04/08/2020 RRTMG https://github.com/NCAR/fv3atm/pull/33 -> https://github.com/NCAR/fv3atm/pull/38
    git cherry-pick 1bae50c9fd5a365226088d861b9d80af9b60de7d
    git cherry-pick ed8347eb0f001de2f53989054c024ad24394a483
    git cherry-pick f1cad8833735ce833bc2ed08a318b77c04e3f536
    git cherry-pick d9841c8d353c262c5aa8e0394ab3ed02542ca8f6
    git cherry-pick 3c64b5a3eaa7608edef1d7c91066d4265e442b27

ufs-weather-model
-----------------

04/06/2020 MONINEDMF https://github.com/NCAR/ufs-weather-model/pull/15 -> https://github.com/NCAR/ufs-weather-model/pull/35
    git cherry-pick e6fe22ffde6425800b7ed2e1cf3748b606806d28
    git cherry-pick a99d5e73f94f5ffb3bd9bbb1a3dc04171026b6d7
    git cherry-pick 5570528f22a9815437e0ef5b413435a8b9fd5881
    git cherry-pick 2bb15be8cc15c69c178b6e443c0ab8c0d75ef318

04/08/2020 RRTMG https://github.com/NCAR/ufs-weather-model/pull/30 -> https://github.com/NCAR/ufs-weather-model/pull/36
    git cherry-pick 5ab446fe423d9a7ba0172ec782e6346c47c98116
    git cherry-pick 1cbc515c6dff00e729523ff6f922e046179d7dcf
    git cherry-pick c6f51f5090411c737a8d79a27bbf08b8d9b00c28
    git cherry-pick 2a43ab170916fbe60c3d918ff16a510928964e51
    git cherry-pick 73f5dcbc93ae297a6a9880988a8f7a7109f99c2e

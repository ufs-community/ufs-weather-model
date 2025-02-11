.. role:: raw-html(raw)
    :format: html

.. _hsd:

********************************************
Hierarchical System Development (HSD) Cases
********************************************

Hierarchical System Development is the ability to engage in development and testing at multiple levels of complexity in numerical weather prediction (NWP) software (such as the :term:`UFS`). It typically includes multiple entry points into development (e.g., atmospheric physics, ocean and ice dynamics, or data assimilation for land models and other earth system components), and it can include both operationally relevant and idealized configurations. 

Although the UFS Weather Model (WM) can be run in any of several configurations, from a single-component atmospheric model to a fully coupled model with multiple earth system components (e.g., atmosphere, ocean, sea-ice, land, mediator), this chapter documents just a few of the cases designed to support hierarchical system development (HSD) within the UFS. Additionally, it explains how to run cases in a container when users do not have access to NOAA :term:`RDHPC systems <RDHPCS>`. 

.. toctree::
   :maxdepth: 3

   CAPE2020
   baroclinic_wave
   HSDcontainer

Currently, users can find information on running the following HSD cases:

   * The :ref:`July 2020 CAPE Case <cape-2020>`
   * The :ref:`Baroclinic Instability Case <baroclinic-wave>`

For a full list of supported WM configurations, view the `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`_ file.

.. attention::

   This chapter is a work in progress. There are a multitude of options for configuring the UFS WM, 
   and this chapter merely details a few supported configurations. It will be expanded over time
   to include a wide variety of idealized test cases for use in research and testing. 

:orphan:

*******************************************************************
Sample ``nems.configure`` File for the ``ATMAQ`` WM Configuration
*******************************************************************

.. code-block:: console

        EARTH_component_list: ATM AQM
        EARTH_attributes::
          Verbosity = 0
        ::
        
        # ATM #
        ATM_model:                      fv3
        ATM_petlist_bounds:             0 271
        ATM_attributes::
          Verbosity = 0
        ::
        
        # AQM #
        AQM_model:                      aqm
        AQM_petlist_bounds:             0 271
        AQM_attributes::
          Verbosity = 0
        ::
        
        # Run Sequence #
        runSeq::
          @180
            ATM phase1
            ATM -> AQM
            AQM
            AQM -> ATM
            ATM phase2
          @
        ::


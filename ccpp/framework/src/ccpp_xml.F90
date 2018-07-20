!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief XML functions and subroutines module.
!!
!! @details The XML module provides functions and 
!!          subroutines for accessing the C versions
!!          of the functions.
!
module ccpp_xml

    use, intrinsic :: iso_c_binding
    use            :: ccpp_types,                                     &
                      only: ccpp_suite_t, ccpp_group_t,               &
                            ccpp_subcycle_t, ccpp_scheme_t,           &
                            CCPP_STR_LEN
    use            :: ccpp_strings,                                   &
                      only: ccpp_fstr, ccpp_cstr, ccpp_free
    use            :: ccpp_errors,                                    &
                      only: ccpp_error, ccpp_warn


    implicit none

    private
    public :: ccpp_xml_load,                                          &
              ccpp_xml_unload,                                        &
              ccpp_xml_ele_find,                                      &
              ccpp_xml_ele_next,                                      &
              ccpp_xml_parse,                                         &
              CCPP_XML_ELE_SUITE,                                     &
              CCPP_XML_ELE_INIT,                                      &
              CCPP_XML_ELE_FINALIZE,                                  &
              CCPP_XML_ELE_GROUP,                                     &
              CCPP_XML_ELE_SUBCYCLE,                                  &
              CCPP_XML_ELE_SCHEME

    interface ccpp_xml_parse
        module procedure ccpp_xml_parse_suite,                        &
                         ccpp_xml_parse_group,                        &
                         ccpp_xml_parse_subcycle,                     &
                         ccpp_xml_parse_fptr
    end interface ccpp_xml_parse

    !>
    !! @brief XML tags for a suite file.
    !!
    !! @details These suite xml tags must match the elements and attributes
    !!          of the suite.xsd.
    !
    character(len=*), parameter :: CCPP_XML_ELE_SUITE    = "suite"
    character(len=*), parameter :: CCPP_XML_ELE_INIT     = "init"
    character(len=*), parameter :: CCPP_XML_ELE_FINALIZE = "finalize"
    character(len=*), parameter :: CCPP_XML_ELE_GROUP    = "group"
    character(len=*), parameter :: CCPP_XML_ELE_SUBCYCLE = "subcycle"
    character(len=*), parameter :: CCPP_XML_ELE_SCHEME   = "scheme"

    character(len=*), parameter :: CCPP_XML_ATT_NAME     = "name"
    character(len=*), parameter :: CCPP_XML_ATT_LOOP     = "loop" 
    character(len=*), parameter :: CCPP_XML_ATT_LIB      = "lib"
    character(len=*), parameter :: CCPP_XML_ATT_VER      = "ver"

    interface
        integer(c_int32_t)                                             &
        function ccpp_xml_load                                         &
                 (filename, xml, root)                                 &
                 bind(c, name='ccpp_xml_load')
            import :: c_int32_t, c_char, c_ptr
            character(kind=c_char), dimension(*) :: filename
            type(c_ptr)                          :: xml
            type(c_ptr)                          :: root
        end function ccpp_xml_load

        integer(c_int32_t)                                             &
        function ccpp_xml_unload                                       &
                 (xml)                                                 &
                 bind(c, name='ccpp_xml_unload')
            import :: c_int32_t, c_ptr
            type(c_ptr)                          :: xml
        end function ccpp_xml_unload

        integer(c_int32_t)                                             &
        function ccpp_xml_ele_find_c                                   &
                 (xml, name, ele)                                      &
                 bind(c, name='ccpp_xml_ele_find')
            import :: c_int32_t, c_ptr, c_char
            type(c_ptr)                          :: xml
            character(kind=c_char), dimension(*) :: name
            type(c_ptr)                          :: ele
        end function ccpp_xml_ele_find_c

        integer(c_int32_t)                                             &
        function ccpp_xml_ele_next_c                                   &
                 (xml, name, ele)                                      &
                 bind(c, name='ccpp_xml_ele_next')
            import :: c_int32_t, c_ptr, c_char
            type(c_ptr)                          :: xml
            character(kind=c_char), dimension(*) :: name
            type(c_ptr)                          :: ele
        end function ccpp_xml_ele_next_c

        integer(c_int32_t)                                             &
        function ccpp_xml_ele_contents                                 &
                 (xml, val)                                            &
                 bind(c, name='ccpp_xml_ele_contents')
            import :: c_int32_t, c_ptr, c_char
            type(c_ptr)                          :: xml
            type(c_ptr)                          :: val
        end function ccpp_xml_ele_contents

        integer(c_int32_t)                                             &
        function ccpp_xml_ele_count                                    &
                 (xml, name, n)                                        &
                 bind(c, name='ccpp_xml_ele_count')
            import :: c_int32_t, c_ptr, c_char
            type(c_ptr)                          :: xml
            character(kind=c_char), dimension(*) :: name
            integer(c_int32_t)                   :: n
        end function ccpp_xml_ele_count

        integer(c_int32_t)                                             &
        function ccpp_xml_ele_att                                      &
                 (node, name, val)                                     &
                 bind(c, name='ccpp_xml_ele_att')
            import :: c_int32_t, c_ptr, c_char
            type(c_ptr)                          :: node
            character(kind=c_char), dimension(*) :: name
            type(c_ptr)                          :: val
        end function ccpp_xml_ele_att

    end interface

    contains

    !>
    !! Find an element in a XML structure.
    !!
    !! @param[in    ] xml      The xml structure.
    !! @param[in,out] name     The element name to find.
    !! @param[   out] ele      The element (if found).
    !! @param[   out] ierr     Integer error flag.
    !
    subroutine ccpp_xml_ele_find(xml, name, ele, ierr)
        type(c_ptr),            intent(in   ) :: xml
        character(len=*),       intent(in   ) :: name
        type(c_ptr),            intent(  out) :: ele
        integer,                intent(  out) :: ierr

        ierr = ccpp_xml_ele_find_c(xml, ccpp_cstr(name), ele)
    end subroutine ccpp_xml_ele_find

    !>
    !! Move to the next occurance of an element in a
    !! XML structure.
    !!
    !! @param[in    ] xml      The xml structure.
    !! @param[in,out] name     The element name to find.
    !! @param[   out] ele      The element (if found).
    !! @param[   out] ierr     Integer error flag.
    !
    subroutine ccpp_xml_ele_next(xml, name, ele, ierr)
        type(c_ptr),            intent(inout) :: xml
        character(len=*),       intent(in   ) :: name
        type(c_ptr),            intent(inout) :: ele
        integer,                intent(  out) :: ierr

        ierr = ccpp_xml_ele_next_c(xml, ccpp_cstr(name), ele)
    end subroutine ccpp_xml_ele_next

    !>
    !! Parse a suite element from an XML structure.
    !!
    !! @param[in    ] node     The current xml node.
    !! @param[in,out] suite    The ccpp_suite_t type to parse into.
    !! @param[   out] ierr     Integer error flag.
    !
    subroutine ccpp_xml_parse_suite(node, suite, ierr)
        type(c_ptr),                intent(in   ) :: node
        type(ccpp_suite_t),         intent(inout) :: suite
        integer,                    intent(  out) :: ierr

        type(c_ptr), target                       :: tmp

        tmp = c_null_ptr

        ! Get the suite name
        ierr = ccpp_xml_ele_att(node, ccpp_cstr(CCPP_XML_ATT_NAME), tmp)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve suite name')
            return
        end if
        suite%name = ccpp_fstr(tmp)
        call ccpp_free(tmp)

        tmp = c_null_ptr

        ! Get the optional library name
        ierr = ccpp_xml_ele_att(node, ccpp_cstr(CCPP_XML_ATT_LIB), tmp)
        if (ierr == 0) then
            suite%library = ccpp_fstr(tmp)
            call ccpp_free(tmp)
            tmp = c_null_ptr
        else
            suite%library = suite%name
        end if

        ! Get the optional library version
        ierr = ccpp_xml_ele_att(node, ccpp_cstr(CCPP_XML_ATT_VER), tmp)
        if (ierr == 0) then
            suite%version = ccpp_fstr(tmp)
            call ccpp_free(tmp)
            tmp = c_null_ptr
        else
            allocate(character(CCPP_STR_LEN) :: suite%version, stat=ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to allocate suite library version')
                return
            end if
            suite%version = ''
            ierr = 0
        end if

        ! Count the number of groups
        ierr = ccpp_xml_ele_count(node, ccpp_cstr(CCPP_XML_ELE_GROUP), suite%groups_max)
        if (ierr /= 0) then
            call ccpp_error('Unable count the number of groups')
            return
        end if

        allocate(suite%groups(suite%groups_max), stat=ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to allocate groups')
            return
        end if

    end subroutine ccpp_xml_parse_suite

    !! Group parsing from an XML file.
    !!
    !! @param[in    ] node       The current xml node.
    !! @param[in    ] max_groups The maximum number of groups.
    !! @param[in,out] group      The ccpp_group_t type to parse into.
    !! @param[   out] ierr       Integer error flag.
    !
    subroutine ccpp_xml_parse_group(node, max_groups, group, ierr)
        type(c_ptr),                intent(in   ) :: node
        integer,                    intent(in   ) :: max_groups
        type(ccpp_group_t),         intent(inout) :: group
        integer,                    intent(  out) :: ierr

        type(c_ptr), target                       :: tmp
        character(kind=c_char), target            :: stmp(CCPP_STR_LEN)

        tmp = c_null_ptr

        ! Get the group name
        ierr = ccpp_xml_ele_att(node, ccpp_cstr(CCPP_XML_ATT_NAME), tmp)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve group name')
            return
        end if
        group%name = ccpp_fstr(tmp)
        call ccpp_free(tmp)

        tmp = c_null_ptr

        ! Count the number of subcycles in this group
        ierr = ccpp_xml_ele_count(node, ccpp_cstr(CCPP_XML_ELE_SUBCYCLE), &
                                  group%subcycles_max)
        if (ierr /= 0) then
            call ccpp_error('Unable to count the number of: ' // &
                            CCPP_XML_ELE_SUBCYCLE)
            return
        end if

        allocate(group%subcycles(group%subcycles_max), stat=ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to allocate subcycles')
            return
        end if

    end subroutine ccpp_xml_parse_group

    !>
    !! Subcycle parsing from an XML file.
    !!
    !! @param[in    ] node          The current xml node.
    !! @param[in    ] max_subcycles The maximum number of subcycles.
    !! @param[in,out] subcycle      The ccpp_subcycle_t type to parse into.
    !! @param[   out] ierr          Integer error flag.
    !
    subroutine ccpp_xml_parse_subcycle(node, max_subcycles, subcycle, ierr)
        type(c_ptr),                intent(in   ) :: node
        integer,                    intent(in   ) :: max_subcycles
        type(ccpp_subcycle_t),      intent(inout) :: subcycle
        integer,                    intent(  out) :: ierr

        type(c_ptr), target                       :: tmp
        character(kind=c_char), target            :: stmp(CCPP_STR_LEN)


        tmp = c_null_ptr

        ! Get the subcycle loop number
        ierr = ccpp_xml_ele_att(node, ccpp_cstr(CCPP_XML_ATT_LOOP), tmp)
        if (ierr /= 0) then
            call ccpp_error('Unable to find subcycle attribute: ' // CCPP_XML_ATT_LOOP)
            return
        else
            stmp = ccpp_fstr(tmp)
            read(stmp, *, iostat=ierr) subcycle%loop
            call ccpp_free(tmp)
            tmp = c_null_ptr
            if (ierr /= 0) then
                call ccpp_error('Unable to convert subcycle attribute "' // &
                                 CCPP_XML_ATT_LOOP // '" to an integer')
            return
            end if
        end if

        ! Count the number of schemes
        ierr = ccpp_xml_ele_count(node, ccpp_cstr(CCPP_XML_ELE_SCHEME), &
                                  subcycle%schemes_max)
        if (ierr /= 0) then
            call ccpp_error('Unable to count the number of: ' // &
                            CCPP_XML_ELE_SCHEME)
            return
        end if

        allocate(subcycle%schemes(subcycle%schemes_max), stat=ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to allocate subcycles')
            return
        end if

    end subroutine ccpp_xml_parse_subcycle

    !>
    !! Function pointer (scheme/init/finalize) parsing from an XML file.
    !!
    !! @param[in    ] node     The current xml node.
    !! @param[in    ] lib      The default library name.
    !! @param[in    ] ver      The default library version.
    !! @param[in,out] fptr     The ccpp_scheme_t type to load into.
    !! @param[   out] ierr     Integer error flag.
    !
    subroutine ccpp_xml_parse_fptr(node, lib, ver, fptr, ierr)
        type(c_ptr),                intent(in   ) :: node
        character(len=*),           intent(in   ) :: lib
        character(len=*),           intent(in   ) :: ver
        type(ccpp_scheme_t),        intent(inout) :: fptr
        integer,                    intent(  out) :: ierr

        type(c_ptr), target                       :: tmp

        tmp = c_null_ptr

        ierr = ccpp_xml_ele_contents(node, tmp)
        if (ierr /= 0) then
            return
        end if

        fptr%name = ccpp_fstr(tmp)
        call ccpp_free(tmp)

        ierr = ccpp_xml_ele_att(node, ccpp_cstr(CCPP_XML_ATT_LIB), tmp)
        if (ierr == 0) then
            fptr%library = ccpp_fstr(tmp)
            call ccpp_free(tmp)
        else
            fptr%library = lib
        end if

        ierr = ccpp_xml_ele_att(node, ccpp_cstr(CCPP_XML_ATT_VER), tmp)
        if (ierr == 0) then
            fptr%version = ccpp_fstr(tmp)
            call ccpp_free(tmp)
        else
            fptr%version = ver
        end if

        ierr = 0
    end subroutine ccpp_xml_parse_fptr

end module ccpp_xml

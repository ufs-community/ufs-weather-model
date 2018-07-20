/*
 * This work (Common Community Physics Package), identified by NOAA, NCAR,
 * CU/CIRES, is free of known copyright restrictions and is placed in the
 * public domain.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file ccpp_xml.h
 *
 * Routines and functions for processing xml files.
 *
 * @ingroup XML
 * @{
 **/
#ifndef CCPP_XML_H
#define CCPP_XML_H

#ifdef __cplusplus
extern "C"
{
#endif

/** Load a XML file. **/
int ccpp_xml_load(const char *, void **, void **);

/** Unload the XML document and finish using the XML library **/
int ccpp_xml_unload(void **);

/** Find the first occurance of the specified element within in a XML
 *  document/node **/
int ccpp_xml_ele_find(void **, const char *, void **);

/** Find the next occurance of the specified element within in a XML node **/
int ccpp_xml_ele_next(void **, const char *, void **);

/** Get the contents of the node **/
int ccpp_xml_ele_contents(void **, char **);

/** Count the number of specifid elements within in a XML document/node **/
int ccpp_xml_ele_count(void **, const char *, int *);

/** Get the attribute at the node. **/
int ccpp_xml_ele_att(void **, const char *, char **);

#ifdef __cplusplus
}                               /* extern "C" */
#endif

#endif                          /* CCPP_XML_H */

/**
 * @}
 **/

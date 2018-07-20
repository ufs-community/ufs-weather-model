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
 * @file ccpp_xml.c
 *
 * @brief Routines and functions for processing a XML file.
 *        This is a very thin layer around libxml2.
 *
 *
 * @ingroup XML
 * @{
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <sysexits.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "ccpp_xml.h"

/**
 * Read a xml file and load the information.
 *
 * @param[in]  filename  The xml file name.
 * @param[out] xml       The xml document pointer.
 * @param[out] root      The root node of the xml document.
 * @retval     0         If it was sucessful.
 * @retval     1         If there was an error.
 **/
int
ccpp_xml_load(const char *filename, void **xml, void **root)
{

	/* Read the file into a document tree */
	*xml = (void *)xmlReadFile(filename, NULL, 0);
	if (*xml == NULL) {
		warnx("Failed to parse %s", filename);
		return(EXIT_FAILURE);
	}

	*root = (void *)xmlDocGetRootElement((xmlDocPtr)(*xml));

	return(EXIT_SUCCESS);
}

/**
 * Unload the XML document and clean-up the XML parser.
 *
 * @param[in] xml   The xml document pointer.
 * @retval     0    If it was sucessful.
 * @retval     1    If there was an error.
 **/
int
ccpp_xml_unload(void **xml)
{
	xmlDocPtr doc = NULL;         /**< XML document tree **/

	doc = (xmlDocPtr)(*xml);

	/* Free the document tree */
	xmlFreeDoc(doc);

	/* Clean up the parser */
	xmlCleanupParser();

	return(EXIT_SUCCESS);
}

/**
 * Get the first occurance of the node.
 *
 * @param[in]  node    The toplevel node pointer to start from.
 * @param[in]  name    The name element to retrieve.
 * @param[out] ele     The first occurance of the element.
 * @retval     0       If it was sucessful.
 * @retval     1       If there was an error.
 **/
int
ccpp_xml_ele_find(void **node, const char *name, void **ele)
{
	xmlNodePtr cur = NULL;         /**< XML tree root node **/

	cur = (xmlNodePtr)(*node);

	/* Loop through all children finding the first requested element*/
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if (xmlStrcmp(cur->name, (const xmlChar *)name) == 0) {
			*ele = (void *)cur;
			break;
		}
		cur = cur->next;
	}
	if (!*ele) {
		return(EXIT_FAILURE);
	}

	return(EXIT_SUCCESS);
}

/**
 * Get the next occurance of the node.
 *
 * This uses xmlNextElementSibling() followed by a check of
 * the name.
 *
 * @param[in]  node    The toplevel node pointer to start from.
 * @param[in]  name    The name element to retrieve.
 * @param[out] ele     The next occurance of the element.
 * @retval     0       If it was sucessful.
 * @retval     1       If there was an error.
 **/
int
ccpp_xml_ele_next(void **node, const char *name, void **ele)
{
	xmlNodePtr cur = NULL;         /**< XML tree root node **/

	cur = (xmlNodePtr)(*node);

	cur = xmlNextElementSibling(cur);
	/* Loop through all siblings finding the element requested */
	while (cur != NULL) {
		if (xmlStrcmp(cur->name, (const xmlChar *)name) == 0) {
			*ele = (void *)cur;
			break;
		}
		cur = xmlNextElementSibling(cur);
	}
	if (!*ele) {
		return(EXIT_FAILURE);
	}

	return(EXIT_SUCCESS);
}

/**
 * Count the number of elements within the XML node.
 *
 * @param[in]  node    The toplevel node pointer to start from.
 * @param[in]  name    The name element to count.
 * @param[out] n       The number of times the element was found.
 * @retval     0       If it was sucessful.
 * @retval     1       If there was an error.
 **/
int
ccpp_xml_ele_count(void **node, const char *name, int *n)
{
	xmlNodePtr cur = NULL;         /**< XML tree root node **/

	cur = (xmlNodePtr)(*node);

	*n = 0;

	/* Count the number of elements */
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if (xmlStrcmp(cur->name, (const xmlChar *)name) == 0) {
			++(*n);
		}
		cur = cur->next;
	}

	return(EXIT_SUCCESS);
}

/**
 * Get the contents of a node.
 *
 * @param[in]  node    The toplevel node pointer to start from.
 * @param[out] value   The value of the attribute.
 * @retval     0       If it was sucessful.
 * @retval     1       If there was an error.
 **/
int
ccpp_xml_ele_contents(void **node, char **value)
{
	int n          = 0;            /**< String length **/
	xmlNodePtr cur = NULL;         /**< XML tree node **/
	xmlChar *tmp   = NULL;         /**< The contents value **/

	cur = (xmlNodePtr)(*node);

	tmp = xmlNodeGetContent(cur);
	if (!tmp) {
		return(EXIT_FAILURE);
	}

	n = strlen((char *)tmp);
	*value = malloc((n+1) * sizeof(char));
	strncpy(*value, (char *)tmp, n * sizeof(char));
	(*value)[n] = '\0';
	xmlFree(tmp);

	return(EXIT_SUCCESS);
}

/**
 * Get the attribute at the node.
 *
 * @param[in]  node    The toplevel node pointer to start from.
 * @param[in]  name    The name of the attribute to get.
 * @param[out] value   The value of the attribute.
 * @retval     0       If it was sucessful.
 * @retval     1       If there was an error.
 **/
int
ccpp_xml_ele_att(void **node, const char *name, char **value)
{
	int n          = 0;            /**< String length **/
	xmlNodePtr cur = NULL;         /**< XML tree node **/
	xmlChar *tmp   = NULL;         /**< The attribute value **/

	cur = (xmlNodePtr)(*node);

	tmp = xmlGetProp(cur, (const xmlChar *)name);
	if (!tmp) {
		return(EXIT_FAILURE);
	}

	n = strlen((char *)tmp);
	*value = malloc((n+1) * sizeof(char));
	strncpy(*value, (char *)tmp, n * sizeof(char));
	(*value)[n] = '\0';
	xmlFree(tmp);

	return(EXIT_SUCCESS);
}

/**
 * @}
 **/

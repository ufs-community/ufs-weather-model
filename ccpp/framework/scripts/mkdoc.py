#!/usr/bin/env python

#
# Functions to generate basic documentation in HTML and LaTeX for CCPP metadata
#

# DH* TODO: create a Python module metadata.py with a class Metadata
# and use this for ccpp_prebuild.py; create to_html and to_latex routines for it

import logging

from common import decode_container, escape_tex

###############################################################################

def metadata_to_html(metadata, model, filename):
    """Create an HTML page with a table that lists each variable provided
    by the model. Contrary to metadata_to_latex below, this table does not
    include information on variables requested by schemes. The primary use
    of the HTML table is to help physics scheme developers to identify the
    variables they need when writing a CCPP-compliant scheme."""

    shading = { 0 : 'darkgray', 1 : 'lightgray' }
    success = True

    # Header
    html = '''<html>
<title>CCPP variables provided by model {model}</title>
<body>
<h1>CCPP variables provided by model {model}</h1>
<table columns="6" border="1" cellpadding="4">
<tr>
  <th bgcolor="{bgcolor}" align="left" >standard_name</th>
  <th bgcolor="{bgcolor}" align="left" >long_name    </th>
  <th bgcolor="{bgcolor}" align="left" >units        </th>
  <th bgcolor="{bgcolor}" align="right">rank         </th>
  <th bgcolor="{bgcolor}" align="left" >type         </th>
  <th bgcolor="{bgcolor}" align="left" >kind         </th>
  <th bgcolor="{bgcolor}" align="left" >source       </th>
  <th bgcolor="{bgcolor}" align="left" >{model} name </th>
</tr>
'''.format(model=model, bgcolor = shading[0])

    count = 0
    for var_name in sorted(metadata.keys()):
        for var in metadata[var_name]:
            # Alternate shading, count is 0 1 0 1 ...
            count = (count+1) % 2
            # ... create html row ...
            line = '''<tr>
  <td bgcolor="{bgcolor}" align="left" >{v.standard_name}</td>
  <td bgcolor="{bgcolor}" align="left" >{v.long_name}    </td>
  <td bgcolor="{bgcolor}" align="left" >{v.units}        </td>
  <td bgcolor="{bgcolor}" align="right">{rank}           </td>
  <td bgcolor="{bgcolor}" align="left" >{v.type}         </td>
  <td bgcolor="{bgcolor}" align="left" >{v.kind}         </td>
  <td bgcolor="{bgcolor}" align="left" >{container}      </td>
  <td bgcolor="{bgcolor}" align="left" >{v.local_name}   </td>
</tr>'''.format(v=var, rank=var.rank.count(':'), container = decode_container(var.container), bgcolor=shading[count])
            html += line

    # Footer
    html += '''</table>
</body>
</html>
'''

    with open(filename, 'w') as f:
        f.write(html)

    logging.info('Metadata table for model {0} written to {1}'.format(model, filename))
    return success


def metadata_to_latex(metadata_define, metadata_request, model, filename):
    """Create a LaTeX document with a table that lists  each variable provided
    and/or requested. Uses the GMTB LaTeX templates and style definitons in gmtb.sty."""

    shading = { 0 : 'darkgray', 1 : 'lightgray' }
    success = True

    var_names = sorted(list(set(metadata_define.keys() + metadata_request.keys())))

    latex = '''\\documentclass[12pt,letterpaper,oneside]{{scrbook}}

\\usepackage{{import}}
\\import{{../common/}}{{gmtb.sty}}
\\renewcommand{{\\thesection}}{{\\arabic{{section}}}}
\\renewcommand{{\\thesubsection}}{{\\arabic{{section}}.\\arabic{{subsection}}}}

\\begin{{document}}

\\section{{CCPP variables provided by model {model} vs requested by pool of physics}}\\label{{sec_ccpp_variables}}
\\subsection{{List of variables}}
\\begin{{longtable}}{{l}}'''.format(model=model)

    for var_name in var_names:
        if var_name in metadata_define.keys():
            var = metadata_define[var_name][0]
        else:
            var = metadata_request[var_name][0]
        line = '''
\hyperlink{{{standard_name_ref}}}{{\\blue\\underline{{\\execout{{{standard_name}}}}}}} \\\\'''.format(
                     standard_name=escape_tex(var.standard_name), standard_name_ref=var.standard_name)
        latex += line

    latex += '''
\\end{longtable}\\pagebreak
\\subsection{Description of variables}
{{\\small\\begin{description}
'''

    for var_name in var_names:
        if var_name in metadata_define.keys():
            var = metadata_define[var_name][0]
            target = escape_tex(decode_container(var.container))
            local_name = escape_tex(var.local_name)
        else:
            var = metadata_request[var_name][0]
            target = 'MISSING'
            local_name = 'MISSING'
        if var_name in metadata_request.keys():
            requested_list = [ escape_tex(decode_container(v.container)) for v in metadata_request[var_name] ]
            # for the purpose of the table, just output the name of the subroutine
            for i in xrange(len(requested_list)):
                entry = requested_list[i]
                requested_list[i] = entry[entry.find('SUBROUTINE')+len('SUBROUTINE')+1:]
            requested = '\\newline '.join(sorted(requested_list))
        else:
            requested = 'NOT REQUESTED'

        # Create output
        text = '''
\\begin{{samepage}}\\item{{
\hypertarget{{{standard_name_ref}}}{{\\blue\\exec{{{standard_name}}}}}}}\\\\ \\nopagebreak
\\begin{{tabular}}{{ll}}
\\execout{{long\_name }} & \\execout{{{long_name}          }} \\\\
\\execout{{units      }} & \\execout{{{units}              }} \\\\
\\execout{{rank       }} & \\execout{{{rank}               }} \\\\
\\execout{{type       }} & \\execout{{{type}               }} \\\\
\\execout{{kind       }} & \\execout{{{kind}               }} \\\\
\\execout{{source     }} & \\execout{{{target}             }} \\\\
\\execout{{local\_name}} & \\execout{{{local_name}         }} \\\\
\\execout{{requested  }} & \\execout{{\\vtop{{{requested}}}}} \\\\
\\end{{tabular}}
\\vspace{{4pt}}
\\end{{samepage}}'''.format(standard_name=escape_tex(var.standard_name), standard_name_ref=var.standard_name,
                          long_name=escape_tex(var.long_name),
                          units=escape_tex(var.units),
                          rank=var.rank.count(':'),
                          type=escape_tex(var.type),
                          kind=escape_tex(var.kind),
                          target=target,
                          local_name=local_name,
                          requested=requested)
        latex += text
    # Footer
    latex += '''
\\end{description}}}
\\end{document}
'''

    with open(filename, 'w') as f:
        f.write(latex)

    logging.info('Metadata table for model {0} written to {1}'.format(model, filename))
    return success

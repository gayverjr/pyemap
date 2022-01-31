## Copyright (c) 2017-2022, James Gayvert, Ruslan Tazhigulov, Ksenia Bravaya
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from .data import res_name_to_char


def extract_resname(residue):
    try:
        resname = residue.resname
        if resname.upper() in res_name_to_char:
            return resname
        resnum = str(residue.full_id[3][1])
        return resname[:resname.rfind(resnum)]
    except Exception:
        return resname


def validate_binary_params(dist_def, edge_prune, sdef):
    ''' Checks if dist_def, sdef, and edge_prune are defined properly. Accepts old 0,1 inputs as well for compatibility with older versions/eMap.
    '''
    try:
        if int(dist_def) == 0:
            dist_def = 'COM'
        elif int(dist_def) == 1:
            dist_def = 'CATM'
        else:
            raise PyeMapGraphException("Improper specification of dist_def. Should be 'COM' or 'CATM'.")
    except ValueError:
        if dist_def.upper() not in ['COM', 'CATM']:
            raise PyeMapGraphException("Improper specification of dist_def. Should be 'COM' or 'CATM'.")
        else:
            dist_def = dist_def.upper()
    try:
        if sdef is None:
            pass
        elif int(sdef) == 0:
            sdef = 'RD'
        elif int(sdef) == 1:
            sdef = 'RSA'
        else:
            sdef = None
    except ValueError:
        if sdef.upper() not in ['RD', 'RSA']:
            sdef = None
        else:
            sdef = sdef.upper()
    try:
        if int(edge_prune) == 0:
            edge_prune = 'DEGREE'
        elif int(edge_prune) == 1:
            edge_prune = 'PERCENT'
        else:
            raise PyeMapGraphException("Improper specification of edge_prune. Should be 'DEGREE' or 'PERCENT'.")
    except ValueError:
        if edge_prune.upper() not in ['DEGREE', 'PERCENT']:
            raise PyeMapGraphException("Improper specification of edge_prune. Should be 'DEGREE' or 'PERCENT'.")
        else:
            edge_prune = edge_prune.upper()
    return dist_def, edge_prune, sdef

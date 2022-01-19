from .data import res_name_to_char

def extract_resname(residue):
    try:
        resname = residue.resname
        if resname.upper() in res_name_to_char:
            return resname
        resnum = str(residue.full_id[3][1])
        return resname[:resname.rfind(resnum)]
    except:
        return resname


def validate_binary_params(dist_def,edge_prune,sdef):
    ''' Checks if dist_def, sdef, and edge_prune are defined properly. Accepts old 0,1 inputs as well for compatibility with older versions/eMap.
    '''
    try:
        if int(dist_def)==0:
            dist_def = 'COM'
        elif int(dist_def) == 1:
            dist_def = 'CATM'
        else:
            raise PyeMapGraphException("Improper specification of dist_def. Should be 'COM' or 'CATM'.")
    except ValueError:
        if dist_def.upper() not in ['COM','CATM']:
            raise PyeMapGraphException("Improper specification of dist_def. Should be 'COM' or 'CATM'.")
        else:
            dist_def = dist_def.upper()
    try:
        if sdef==None:
            pass
        elif int(sdef)==0:
            sdef = 'RD'
        elif int(sdef) == 1:
            sdef = 'RSA'
        else:
            sdef = None
    except ValueError:
        if sdef.upper() not in ['RD','RSA']:
            sdef = None
        else:
            sdef = sdef.upper()
    try:
        if int(edge_prune)==0:
            edge_prune = 'DEGREE'
        elif int(edge_prune) == 1:
            edge_prune = 'PERCENT'
        else:
            raise PyeMapGraphException("Improper specification of edge_prune. Should be 'DEGREE' or 'PERCENT'.")
    except ValueError:
        if edge_prune.upper() not in ['DEGREE','PERCENT']:
            raise PyeMapGraphException("Improper specification of edge_prune. Should be 'DEGREE' or 'PERCENT'.")
        else:
            edge_prune = edge_prune.upper()
    return dist_def,edge_prune,sdef
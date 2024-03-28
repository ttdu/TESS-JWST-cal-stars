import pyvo as vo

def get_name_type(ticid):
    service = vo.dal.TAPService("https://simbad.cds.unistra.fr/simbad/sim-tap/")

    text = f"""
    SELECT main_id,
           otype,
           sp_type
    FROM basic JOIN ident ON oidref = oid
    WHERE id = 'TIC {ticid}';
    """
    resultset = service.search(text)
    try:
        name = resultset['main_id'].data[0]
        otype = resultset['otype'].data[0]
        stype = resultset['sp_type'].data[0]
    except:
        name = otype = stype = "???"
    return name, otype, stype

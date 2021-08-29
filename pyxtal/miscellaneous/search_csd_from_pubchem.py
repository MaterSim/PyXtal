import pubchempy as pcp
import urllib, json
from pprint import pprint

def get_similar_cids(base, MaxRecords):
    '''''''''
    Parameters:
    base: PubChem CID of Starting chemical
    MaxRecords: Number of Similar Compounds

    Returns:
    List of the CIDs of PubChem compounds similar to the base compound.

    Accuracy decreases deeper into PubChem search algorithm
    '''''''''
    if type(base) == int: base = str(base)
    cids = pcp.get_compounds(base, searchtype="similarity", MaxRecords=MaxRecords)
    results = []
    for x in cids:
        csd_codes = check_for_ccdc_structures(x.cid)
        if len(csd_codes)>0:
            d = {"cid": x.cid,
                 "smiles": x.canonical_smiles,
                 "name": x.iupac_name,
                 "csd_codes": csd_codes}
            results.append(d)
            pprint(d)

    return results

def check_for_ccdc_structures(cid):
    ''''''''''
    Parameters:
    cid: PubChem cid

    Returns:
    CIDs that have CCDC crystal structure data

    '''''''''
    url0 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'
    cid=str(cid)
    url = url0 + cid + '/JSON'
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())

    csd_codes = []
    
    if len(data['Record']['Section'][0]['Section']) == 3:
        infos = data['Record']['Section'][0]['Section'][2]['Section'][0]['Information']
        for info in infos:
            csd_codes.append(info['Value']['StringWithMarkup'][0]['String'])
    return csd_codes


if __name__ == "__main__":
    s = get_similar_cids(2244, 25) 

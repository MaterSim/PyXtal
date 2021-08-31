from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import pubchempy as pcp
import sys

def get_similar_cids(base,MaxRecords):
    '''''''''
    Parameters:
    base: PubChem CID of Starting chemical
    MaxRecords: Number of Similar Compounds

    Returns:
    List of the CIDs of PubChem compounds similar to the base compound.

    Accuracy decreases deeper into PubChem search algorithm
    '''''''''

    base=str(base)
    cids=pcp.get_compounds(base ,searchtype="similarity",MaxRecords=MaxRecords)
    cids=[x.to_dict()['cid'] for x in cids]
    return cids



def check_for_ccdc_structures(cids):
    ''''''''''
    Parameters:
    cids: list of PubChem cids

    Returns:
    List of the the given CIDs that have CCDC crystal structure data
    '''''''''
    good_list=[]

    for cid in cids:
        cid=str(cid)
        start="https://pubchem.ncbi.nlm.nih.gov/"
        driver=webdriver.chrome.webdriver.WebDriver("/Users/kevinparrish/github/PyXtal/pyxtal/miscellaneous/web_scalping_script/chromedriver.exe")
        driver.get(start)
        time.sleep(2)
        elem=driver.find_element_by_css_selector('input')
        elem.send_keys(cid)
        elem.send_keys(Keys.RETURN)

        time.sleep(3)

        elem=driver.find_element_by_link_text(cid)
        elem.click()

        time.sleep(4)

        page=driver.page_source
        time.sleep(1)
        if "CCDC" in page:
            good_list.append(cid)
        driver.quit()
    return good_list




def ccdcid_scalper(cids):
    ''''''''''
    Parameters:
    cids: list of PubChem cids

    Returns:
    List of lists. Primary axis same dimension as input list of cids.
    Each element is a list of all the CCDC Numbers provided in the PubChem Page for each PubChem CID.
    If there is no crystal structure data on the PubChem Page for a particular CID, there will be an empty list.



    On CCDC website, bulk cif downloads can be done with using a comma separated list of CCDC numbers.
    Flatten the returned list of this function to get that.
    '''''''''
    CCDC_numbers=[]

    for cid in cids:
        cid=str(cid)
        start="https://pubchem.ncbi.nlm.nih.gov/"
        driver_location=sys.argv[1]
        driver=webdriver.chrome.webdriver.WebDriver(executable_path=str(driver_location))
        driver.get(start)
        time.sleep(2)
        elem=driver.find_element_by_css_selector('input')
        elem.send_keys(cid)
        elem.send_keys(Keys.RETURN)

        time.sleep(3)

        elem=driver.find_element_by_link_text(cid)
        elem.click()

        time.sleep(4)

        elems=driver.find_elements_by_partial_link_text("Ccdcid=")
        urls=[elem.text for elem in elems]
        ccdcids=[url[url.index("=")+1:] for url in urls]
        ccdcids=[int(ccdcid) for ccdcid in ccdcids]

        CCDC_numbers.append(ccdcids)
        driver.quit()
    return CCDC_numbers



if __name__ == "__main__":
    s=get_similar_cids(2244,150) #2244 is aspirins PubChem CID
    print('List of CIDs of Similar PubChem Compounds to Aspirin: ,',s)

    print('\nCCDC Numbers listed on PubChem Page for each PubChem CID')
    ccdc_numbers=ccdcid_scalper(s)

    for i,x in enumerate(ccdc_numbers):
        print(s[i],': ',x)

    input=[]
    [input.extend(x) for x in ccdc_numbers]
    print('\nFlattened array for CCDC bulk search and download: ', input)

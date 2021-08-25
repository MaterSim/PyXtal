from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import pubchempy as pcp

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

    This assumes your system has google chrome v.92 installed. If not, change the chromedriver.exe in the folder.
    Drivers can be found on the selenium project website.

    '''''''''
    good_list=[]

    for cid in cids:
        cid=str(cid)
        start="https://pubchem.ncbi.nlm.nih.gov/"
        driver=webdriver.chrome.webdriver.WebDriver()
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


if __name__ == "__main__":
    s=get_similar_cids(2244,5) #2244 is aspirins PubChem CID. First 5 similar compounds are gathered and checked.
    print(s)

    crystal_cids=check_for_ccdc_structures(s)
    print(crystal_cids)

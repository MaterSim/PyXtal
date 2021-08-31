# Automated Crystal Structure Data Scalper
Uses the collective database of PubChem and the Cambridge Crystallographic Data Centre to download cif files

## Set-Up
In order to get the online scalper set up, first ensure that PyXtal is already installed in your environment. Then, execute these shell commands
```bash
cd PyXtal/pyxtal/miscellaneous/web_scalping_script
pip install -r requirements.txt
```

This script relies on the python package [Selenium](https://www.selenium.dev/), a powerful software that gives pythonic power over internet browsing. 
For it to work, the latest version of a WebDriver needs to be made available to the script. 
This script is written to use Google Chrome. To acquire the webdriver, first have google chrome and installed and find the version number.
Then find the corresponding version number in this link here, [https://chromedriver.chromium.org/downloads](https://chromedriver.chromium.org/downloads) and download the appropriate driver for your operating system.
Unzip the executable file, and find the absolute path to the chromedriver.exe file. The script takes the absolute path as an argument

```bash
python3 pubchem_import.py "/Users/kevinparrish/Downloads/chromedriver"
```

### NOTE
Your system may block the the executable. Be sure to allow the executable permission to run.

## Workflow
The script relies on two functionalities: 

`get_similiar_cids` accepts a PubChem Compound CID and a Maximium threshold as arguments. It then uses the PubChem python API `pubchempy` to search through pubchem database the the threshold number of similar structure compounds. It returns a list of the similar compound CIDs.

`ccdcid_scalper` then takes the list and uses selenium to systematically search the PubChem website for these CIDs, where it will identify any potential links to the CCDC website attached to the page. It will grab all the CDC Numbers available, and return a final list of comma separated CDC Numbers that represent CIF files of crystals similar in structure to the starting input.

Then, the comma separated list can be taken and put into the Identifiers field here [https://www.ccdc.cam.ac.uk/structures/Home/EditSearchForm?ccdc-check=52d8db2b85f6380110d1d48a537426b2](https://www.ccdc.cam.ac.uk/structures/Home/EditSearchForm?ccdc-check=52d8db2b85f6380110d1d48a537426b2). Check only the CSD database to search. The resulting list of cif files can be downloaded all at once through the button in the top left.


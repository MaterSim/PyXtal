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
python pubchem_import.py "/Users/kevinparrish/Downloads/chromedriver"
```

### NOTE
Your system may block the the executable. Be sure to allow the executable permission to run.


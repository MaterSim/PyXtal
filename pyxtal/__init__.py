from pyxtal.version import __version__

name = "pyxtal"


def print_logo():
    """
    print out the logo and version information
    """

    print(
        """
             ______       _    _          _   
            (_____ \     \ \  / /        | |   
             _____) )   _ \ \/ / |_  ____| |  
            |  ____/ | | | )  (|  _)/ _  | | 
            | |    | |_| |/ /\ \ |_( (_| | |___
            |_|     \__  /_/  \_\___)__|_|_____)
                   (____/      """
    )
    print("\n")
    print("----------------------(version", __version__, ")----------------------\n")
    print("A Python package for random crystal generation")
    print("The source code is available at https://github.com/qzhu2017/pyxtal")
    print("Developed by Zhu's group at University of Nevada Las Vegas\n\n")


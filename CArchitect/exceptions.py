class Directory_Not_Found(OSError):
    """My defined exception"""
    def __init__(self, directory):
        self.directory = directory
    def __str__(self):
        return "Directory %s does not exist. Please select a valid directory." %(self.directory)


class No_Input_Directory_Provided(ValueError):
    """My defined exception"""
    pass


class No_PDB_Found(ValueError):
    """My defined exception"""
    pass

class Incorrect_Number_Chains(ValueError):
    """My defined exception"""
    pass

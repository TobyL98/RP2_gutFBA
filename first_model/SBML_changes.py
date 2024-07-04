######################
# SBML_changes.py
######################

'''This code modfies two member ID names in the AGORA western diet models
SBML files that commonly did not match the names in the reactions. The models
should work after correcting for these changes '''

import xml.etree.ElementTree as ET
from pathlib import Path

folder_Path = Path('Agora_HighFiber/sbml')

for file in folder_Path.glob('*.xml'):

    tree = ET.parse(file)
    root = tree.getroot()

    # fins the correct XML tag 'member'
    for member in root.iter('{http://www.sbml.org/sbml/level3/version1/groups/version1}member'):
        # pull out the ID of the specific member
        member_ID = member.attrib['{http://www.sbml.org/sbml/level3/version1/groups/version1}idRef']


        # sets the new IDs
        if member_ID == "R_EX_sbt__45__d__40__e__41__":
            new_member_ID = "R_EX_sbt_D__40__e__41__"
            member.set('{http://www.sbml.org/sbml/level3/version1/groups/version1}idRef', new_member_ID)

        if member_ID =="R_EX_glc__40__e__41__":
            new_member_ID = "R_EX_glc_D__40__e__41__"
            member.set('{http://www.sbml.org/sbml/level3/version1/groups/version1}idRef', new_member_ID)

    # outputs to an updated folder
    output_path = Path('Agora_HighFiber_updated') / file.name
    tree.write(output_path)



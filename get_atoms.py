from open_babel import smiles2ConnectionTable




def identify_last_character (table):
    lines = table.split("\n")
    atoms = []
    for line in lines:
        atoms.append(line[-1])
    return atoms

def find_atoms (character_list):
    numbers = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    atoms = []
    for character in character_list:
        if character not in numbers:
            atoms.append(character)
    return atoms

def get_index_hetero_atoms (atoms_list):
    hetero_atoms = []
    for atom in atoms_list:
        if atom is not "C":
            hetero_atoms.append(index(atom))
    return hetero_atoms


if __name__=='__main__':

    _, MDL = smiles2ConnectionTable("O[C@@H](Cc1ccccc1)C(O)=O")

    print(MDL)
    find_atoms(identify_last_character(MDL))
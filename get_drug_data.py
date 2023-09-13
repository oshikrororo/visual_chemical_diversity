import pymysql
# if you use config file, uncomment following line
# from config import user, password, host, port
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_scaffold(smiles: str) -> (str, str):
    mol = Chem.MolFromSmiles(smiles)
    non_generic = MurckoScaffold.GetScaffoldForMol(mol)
    generic = MurckoScaffold.GetScaffoldForMol(MurckoScaffold.MakeScaffoldGeneric(non_generic))
    non_generic = Chem.MolToSmiles(non_generic)
    generic = Chem.MolToSmiles(generic)
    return generic, non_generic


def get_data(cursor) -> dict:
    data = {}
    filter_phase_4 = "SELECT MOLREGNO, MESH_HEADING, EFO_TERM " \
                     "FROM DRUG_INDICATION " \
                     "WHERE MAX_PHASE_FOR_IND = 4;"
    cursor.execute(filter_phase_4)
    rows = cursor.fetchall()

    generic_ids = {}
    generic_count = 0
    non_generic_ids = {}
    non_generic_count = 0

    for indication in rows:
        molregno, mesh, efo = indication

        get_parent = f"SELECT PARENT_MOLREGNO " \
                     f"FROM MOLECULE_HIERARCHY " \
                     f"WHERE MOLREGNO = {molregno};"
        cursor.execute(get_parent)
        molregno, = cursor.fetchone()

        get_id = f"SELECT CHEMBL_ID, PREF_NAME " \
                 f"FROM MOLECULE_DICTIONARY " \
                 f"WHERE MOLREGNO = {molregno};"
        cursor.execute(get_id)
        record_id, name = cursor.fetchone()

        if record_id in data:
            if mesh is not None and 'mesh_indications' in data[record_id]:
                data[record_id]['mesh_indications'] += f'|{mesh}'
            if efo is not None and 'efo_indications' in data[record_id]:
                data[record_id]['efo_indications'] += f'|{efo}'
            continue

        if_organic_drug = f"SELECT ORAL, PARENTERAL " \
                          f"FROM MOLECULE_DICTIONARY " \
                          f"WHERE MOLREGNO = {molregno} AND " \
                          f"(ORAL = 1 OR PARENTERAL = 1) " \
                          f"AND (INORGANIC_FLAG = 0) AND (PRODRUG = 0);"
        cursor.execute(if_organic_drug)
        if not len(cursor.fetchall()):
            continue

        get_smiles = f"SELECT CANONICAL_SMILES " \
                     f"FROM COMPOUND_STRUCTURES " \
                     f"WHERE MOLREGNO = {molregno};"
        cursor.execute(get_smiles)
        smiles = cursor.fetchone()
        if smiles is None:
            continue
        smiles, = smiles

        generic_smiles, non_generic_smiles = get_scaffold(smiles)
        if generic_smiles in generic_ids:
            generic_id = generic_ids[generic_smiles]
        else:
            generic_count += 1
            generic_id = generic_count
            generic_ids[generic_smiles] = generic_id
        if non_generic_smiles in non_generic_ids:
            non_generic_id = non_generic_ids[non_generic_smiles]
        else:
            non_generic_count += 1
            non_generic_id = non_generic_count
            non_generic_ids[non_generic_smiles] = non_generic_id

        data[record_id] = {}
        data[record_id]['name'] = name
        data[record_id]['smiles'] = smiles
        data[record_id]['generic_scaffold'] = generic_smiles
        data[record_id]['generic_id'] = generic_id
        data[record_id]['non_generic_scaffold'] = non_generic_smiles
        data[record_id]['non_generic_id'] = non_generic_id
        if efo is not None:
            data[record_id]['efo_indications'] = efo
        if mesh is not None:
            data[record_id]['mesh_indications'] = mesh

    print(f'Generic scaffolds amount: {generic_count}')
    print(f'Non-generic scaffolds amount: {non_generic_count}')
    return data


def write_tsv(path, data):
    with open(path, 'w') as file:
        file.write('Parent Compound\t'
                   'Preferred Name\t'
                   'SMILES\t'
                   'Non-generic Scaffold SMILES\t'
                   'Non-generic Scaffold ID\t'
                   'Generic Scaffold SMILES\t'
                   'Generic Scaffold ID\t'
                   'MESH Indications\t'
                   'EFO Indications\n')
        for record_id, row in data.items():
            try:
                efo = sorted(list(set(row['efo_indications'].split('|'))), key=lambda x: x.lower())
            except KeyError:
                efo = '-'
            try:
                mesh = sorted(list(set(row['mesh_indications'].split('|'))), key=lambda x: x.lower())
            except KeyError:
                mesh = '-'
            file.write(f'{record_id}\t'
                       f'{row["name"]}\t'
                       f'{row["smiles"]}\t'
                       f'{row["non_generic_scaffold"]}\t'
                       f'SCAFFA{row["non_generic_id"]}\t'
                       f'{row["generic_scaffold"]}\t'
                       f'SCAFFB{row["generic_id"]}\t'
                       f'{"|".join(mesh)}\t'
                       f'{"|".join(efo)}\n')


def main():
    # insert mySQL server parameters here
    connection = pymysql.connect(
        host=host,
        user=user,
        database='chembl_33',
        password=password,
        port=port
    )
    try:
        with connection.cursor() as cursor:
            data = get_data(cursor)
    except Exception as ex:
        print(ex)
        return
    finally:
        connection.close()

    # insert path to table here
    path = 'table.tsv'
    write_tsv(path, data)


if __name__ == '__main__':
    main()

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

import rdkit.Chem as Chem
import rdkit.Chem.Draw

import json

app = FastAPI()

origins = [
    "http://localhost:3000",
    "localhost:3000"
]
# If using VSCode + windows, try using your IP 
# instead (see frontent terminal)
#origins = [
#    "http://X.X.X.X:3000",
#    "X.X.X.X:3000"
#]


app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)


def make_routes(index):
    # TODO: use this method to return routes as a tree data structure.
    # routes are found in the routes.json file

    class Tree():
        def __init__(self,root):
            self.root = root
            self.name = root
            self.children = []
        def createNode(self,obj):
            self.children.append(obj)


    class Node():
        def __init__(self, data):
            self.name = data['name']
            self.attributes = { 
                'catalog_entries_count': len(data['catalog_entries']),
                'reaction_name': data['reaction_name'],
                'smarts_template': data['smarts_template'], 
                }
            self.children = []
        def createNode(self,obj):
            self.children.append(obj)

    
    routesjson = open('app/routes.json')
    routesdata = json.load(routesjson)
    singleroute = routesdata[index]

    Route = Tree('O=C(Cn1nnc2ccccc21)N(Cc1ccsc1)c1ccc(Cl)cc1')

    def create_source_nodes(node, data):
        for reaction in data['reactions']:
            if reaction['target'] == node.name:
                for source in reaction['sources']:
                    sourcedata = next(molecule for molecule in data['molecules'] if molecule['smiles'] == source)
                    newnode = {
                        'name': source,
                        'catalog_entries': sourcedata['catalog_entries'],
                        'reaction_name': reaction['name'],
                        'smarts_template': reaction['smartsTemplate'],
                    }
                    SourceNode = Node(newnode)
                    if not sourcedata['is_building_block']:
                        create_source_nodes(SourceNode, data)
                    node.createNode(SourceNode.__dict__)

    
    create_source_nodes(Route, singleroute)

    singleroute['tree'] = Route.__dict__

    return singleroute


def get_molecule_details(smiles, index):
    routesjson = open('app/routes.json')
    routesdata = json.load(routesjson)
    singleroute = routesdata[index]

    molecule_details = next(molecule for molecule in singleroute['molecules'] if molecule['smiles'] == smiles)
    
    return molecule_details
        


def draw_molecule(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    img = Chem.Draw.MolsToGridImage([mol], molsPerRow=1, useSVG=True)
    return img

@app.get("/", tags=["root"])
async def read_root() -> dict:
    return {
        "message": "Welcome to your app.",
    }

# This endpoint takes a smiles string and returns an svg for that molecule. The svg is returned still in xml formatting, the xml-to-jsx conversion to make it usable by React takes place on the frontend.
@app.get("/molecule", tags=["molecule"])
async def get_molecule(smiles: str) -> dict:
    molecule = draw_molecule(smiles)

    # TODO: return svg image
    moleculeSVG = molecule.replace('\n', ' ').replace('\"', "'").replace('<!-- END OF HEADER -->', '').replace('xmlns:xlink', 'xmlnsXlink').replace('xmlns:rdkit="http://www.rdkit.org/xml"', "").replace("xmlns:rdkit='http://www.rdkit.org/xml'", '').replace('xml:space', "fill='#fff' stroke='#000' x='-50' y='-20' xmlSpace").replace("width='200px' height='200px'", "width='100px' height='100px'").split('?> ')

    return {
        "data": moleculeSVG[1],
    }

# This endpoint takes in an integer that represents the index of a specific route in the routes.json array, and returns the nested tree of the route's component molecules, as well as some attributes like the reaction name, smarts template, and count of catalog entries for each molecule.
@app.get("/routes", tags=["routes"])
async def get_routes(route: int) -> dict:
    routes = make_routes(route)

    return {
        "data": routes,
    }

# This endpoint returns a simplified list of all of the routes in the routes.json array.
# It includes some additional information like the score and the count of all building block molecules.
@app.get("/allroutes")
async def get_routes_list() -> dict:
    routesjson = open('app/routes.json')
    routesdata = json.load(routesjson)
    routes_list = []

    for index, route in enumerate(routesdata):
        route_data = {
            'id': index,
            'score': route['score'],
            'building_blocks': sum(molecule['is_building_block'] == True for molecule in route['molecules'])
        }
        routes_list.append(route_data)

    return {
        "routesList": routes_list
    }

# This endpoint takes a Smiles string and the route index (in reference to the original array in routes.json) and returns the details (catalog_entries and is_building_block) from the molecules array in that route
@app.get("/moleculedetails", tags=["molecule"])
async def get_molecule(smiles: str, index: int) -> dict:
    molecule_details = get_molecule_details(smiles, index)

    return {
        "molecule_details": molecule_details
    }
from pydantic import BaseModel


class Molecule(BaseModel):
    identifier: str
    smiles: str

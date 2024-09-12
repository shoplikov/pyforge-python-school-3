from pydantic import BaseModel

class Molecule(BaseModel):
    identifier: str
    smiles: str

    class Config:    
        orm_mode = True
import PyCrystalField as cef
import numpy as np

def main():
    # pcs=[
    #     [0,0,+4],
    #     [0,0,-4],
    #     [3.46,0,0],
    #     [-3.46,0,0],
    #     [1.73,-3,0],
    #     [-1.73,3,0],
    #     [1.73,3,0],
    #     [-1.73,-3,0]
    # ]
    pcs=[[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]]
    pcs=[
        [-0.01741,+0.31515,+0.08198],
        [+0.33257,+0.01741,+0.08198],
        [-0.31515,-0.33257,+0.08198],
        [+0.33410,+0.31592,-0.08199],
        [-0.31592,+0.01818,-0.08199],
        [-0.01818,-0.33410,-0.08199],
    ]
    lfield=cef.Ligands("Yb3+",pcs)
    pcm=lfield.PointChargeModel(LigandCharge=-1,IonCharge=+3)
    E,V=np.linalg.eig(pcm.H)
    E=np.real(E)-min(np.real(E))
    print()
    for e in np.sort(E):
        print(e)
    return None

if __name__=="__main__":
    main()
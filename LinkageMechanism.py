from sympy import *
from IPython.display import display, Math, Markdown, Latex
init_printing(use_latex='mathjax')
from sympy import diff as D

a,b,c,d,e,f,g,h,x,y,x,z,u,w,k,t,r,l,ℓ = symbols('a b c d e f g h x y x z u w k t r \\ell \\ell')
α,β,γ,δ,θ,φ,φ1,φ2 = symbols('alpha beta gamma delta theta varphi varphi_1 varphi_2')
kα,kβ,kγ,kδ,kθ,kφ = symbols('k_{\\alpha} k_{\\beta} k_{\\gamma} k_{\\delta} k_{\\theta} k_{\\varphi}')
kx,ky,kz,kφ1,kφ2 = symbols('k_x k_y k_z k_{\\varphi_1} k_{\\varphi_2}')
Ck = {α:kα,β:kβ,γ:kγ,δ:kδ,θ:kθ,φ:kφ,x:kx,y:ky,z:kz,φ1:kφ1,φ2:kφ2}
ℓα,ℓβ,ℓγ,ℓδ,ℓθ,ℓφ = symbols('\ell_{\\alpha} \ell_{\\beta} \ell_{\\gamma} \ell_{\\delta} \ell_{\\theta} \ell_{\\varphi}')
ℓx,ℓy,ℓz,ℓφ1,ℓφ2 = symbols('\ell_x \ell_y \ell_z \ell_{\\varphi_1} \ell_{\\varphi_2}')
Ls = {α:ℓα,β:ℓβ,γ:ℓγ,δ:ℓδ,θ:ℓθ,φ:ℓφ,x:ℓx,y:ℓy,z:ℓz,φ1:ℓφ1,φ2:ℓφ2}

θ1,θ2,θ3 = symbols('theta_1 theta_2 theta_3')
x1,x2,y1,y2,z1,z2 = symbols('x_1 x_2 y_1 y_2 z_1 z_2')
αt,βt,γt,δt,θt,φt = symbols('\\dot{\\alpha} \\dot{\\beta} \\dot{\\gamma} \\dot{\\delta} \\dot{\\theta} \\dot{\\varphi}')
θ1t,θ2t,θ3t = symbols('\\dot{\\theta_1} \\dot{\\theta_2} \\dot{\\theta_3}')
x1t,x2t,y1t,y2t,z1t,z2t = symbols('\\dot{x_1} \\dot{x_2} \\dot{y_1} \\dot{y_2} \\dot{z_1} \\dot{z_2}')
xt,yt,zt = symbols('\\dot{x} \\dot{y} \\dot{z}')
Ve = {α:αt,β:βt,γ:γt,δ:δt,θ:θt,θ1:θ1t,θ2:θ2t,θ3:θ3t,φ:φt,x:xt,y:yt,z:zt,x1:x1t,x2:x2t,y1:y1t,y2:y2t,z1:z1t,z2:z2t}

def MecSolve(Cg,Eq):    # Cg,Vg Coords e Veloc Generalizadas
    f = len(Cg)-len(Eq)
    Eq = Matrix(Eq)
    J = Eq.jacobian(Cg[f:])     # Jacobiano do sistema
    F = Eq.jacobian([Cg[:f]])   # Obtenção da matriz F
    K = simplify(-(J**-1)*F)    # Obtenção da matriz K
    if Cg[:f] != 1:
        Ks = Matrix([ Ck[i] for i in Cg[1:] ])
        L = simplify( K.jacobian([Cg[0]]) + K.jacobian([Cg[f:]])*Ks )
    else:
        L = Matrix([])
        Tg = [ i**t for i in Cg ]
        Tg = [ Ve[i] for i in Cg ]
        for i in range(f):          # Obtenção da matriz L
            Lprov = K.col(i).jacobian(Cg[:f])*((1/Tg[i])*Matrix(Tg[:f])) + K.col(i).jacobian(Cg[f:])*((1/Tg[i])*Matrix(Tg[f:]))
            L = L.row_join(Lprov)

    if Cg[:f] != 1:
    # Impresão das matrizes Iniciais
        display( Markdown("### Matrizes F e J") )
        MF = [latex(i) for i in F]
        MJ = [latex(i) for i in J]
        VF = "F = \\begin{Bmatrix}"
        for i in range(len(MF)):
            VF += MF[i]
            if i != len(MF)-1:
                VF += '\\\ '
        VF += "\end{Bmatrix}"
        VJ = "\qquad J = \\begin{bmatrix}"
        tam_MJ = len(MJ); sq_MJ = sqrt(tam_MJ)
        for i in range(0,tam_MJ,4):
            for j in range(sq_MJ):
                VJ += MJ[i+j]
                if j != sq_MJ-1:
                    VJ += ' & '
                j += 1
            if i != tam_MJ-4:
                VJ += '\\\ '
        VJ += "\end{bmatrix}"
        LN = VF + VJ
        display(Markdown('$ %s $' %LN))

    # Impresão dos Coeficientes de Velocidade
        display( Markdown("### Coeficientes de Velocidade") )
        MK = [latex(i) for i in K]
        KS = [latex(i) for i in Ks]
        LN = ""
        for i in range(len(KS)):
            LN += KS[i]+" = "+MK[i]+" \qquad "
        display(Markdown('$ %s $' %LN))

        # Impresão dos Coeficientes da Aceleração
        display( Markdown("### Coeficientes da Aceleração") )
        Ms = Matrix([ Ls[i] for i in Cg[1:] ])
        ML = [latex(i) for i in L]
        MS = [latex(i) for i in Ms]
        LN = ""
        for i in range(len(KS)):
            LN += MS[i]+" = "+ML[i]+" \qquad "
        display(Markdown('$ %s $' %LN))
        return

    return(F, J, K, L)
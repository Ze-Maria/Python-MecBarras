a,b,c,d,e,f,g,h,x,y,x,z,u,w,k,t,r,l,ℓ = symbols('a b c d e f g h x y x z u w k t r \\ell \\ell')
α,β,γ,δ = Function ('alpha')(t), Function ('beta')(t),Function ('gamma')(t),Function ('delta')(t)
θ,φ,x,y,z = Function ('theta')(t), Function ('varphi')(t),Function ('x')(t),Function ('y')(t),Function ('z')(t)
αt,βt,γt,δt,θt,φt = symbols('\\dot{\\alpha} \\dot{\\beta} \\dot{\\gamma} \\dot{\\delta} \\dot{\\theta} \\dot{\\varphi}')
xt,yt,zt = symbols('\\dot{x} \\dot{y} \\dot{z}')
αtt,βtt,γtt = symbols('\\ddot{\\alpha} \\ddot{\\beta} \\ddot{\\gamma}')
δtt,θtt,φtt = symbols('\\ddot{\\delta} \\ddot{\\theta} \\ddot{\\varphi}')
xtt,ytt,ztt = symbols('\\ddot{x} \\ddot{y} \\ddot{z}')
Ve = {α:αt,β:βt,γ:γt,δ:δt,θ:θt,θ1:θ1t,θ2:θ2t,θ3:θ3t,φ:φt,x:xt,y:yt,z:zt,x1:x1t,x2:x2t,y1:y1t,y2:y2t,z1:z1t,z2:z2t}
Ae = {α:αtt,β:βtt,γ:γtt,δ:δtt,θ:θtt,φ:φtt,x:xtt,y:ytt,z:ztt}

def SetupGlobalVars(u,v,sᵢ,xₒ,yₒ):
    global M,U,Xₒ,Ẋₒ,Ẍₒ
    M  = Matrix([[cos(sᵢ), -sin(sᵢ)],[sin(sᵢ), cos(sᵢ)]])
    U  = Matrix([u,v])
    Xₒ = Matrix([xₒ,yₒ])
    Ẋₒ = Matrix([D(xₒ,t),D(yₒ,t)])
    Ẍₒ = Matrix([D(xₒ,t,2),D(yₒ,t,2)])

def P():
    return Xₒ + M*U

def Ṗ(sᵢ):
    Ω = Matrix([[0, -1],[1,0]])
    return Ẋₒ + D(sᵢ,t)*M*Ω*U

def P̈(sᵢ):
    Γ = Matrix([[D(sᵢ,t)**2, D(sᵢ,t,2)],[-D(sᵢ,t,2), D(sᵢ,t)**2]])
    return Ẍₒ - M*Γ*U

def SolveCouplerPoint(u,v,sᵢ,xₒ,yₒ):
    SetupGlobalVars(u,v,sᵢ,xₒ,yₒ)
    Pt = Ṗ(sᵢ); Ptt = P̈(sᵢ)
    for i in Ac:
        Ptt = Ptt.subs( D(i,t,2),Ae[i] )
        Ptt = Ptt.subs( D(i,t),Ve[i] )
    for i in Vc:
        Pt = Pt.subs( D(i,t),Ve[i] )

    desloc = [latex(i) for i in P()]
    veloc  = [latex(i) for i in Pt]
    acel   = [latex(i) for i in Ptt]
    d0 = desloc[0]; d1 = desloc[1]
    v0 = veloc[0]; v1 = veloc[1]
    a0 = acel[0]; a1 = acel[1]
    d0 = d0.replace(" ", ""); d1 = d1.replace(" ", "")
    v0 = v0.replace(" ", ""); v1 = v1.replace(" ", "")
    a0 = a0.replace(" ", ""); a1 = a1.replace(" ", "")
    A = B = C = E = G = H = 1
    while (A > 0 or B > 0 or C > 0 or E > 0 or G > 0 or H > 0):
        A = d0.find('{\\left(t'); B = d1.find('{\\left(t')
        C = v0.find('{\\left(t'); E = v1.find('{\\left(t')
        G = a0.find('{\\left(t'); H = a1.find('{\\left(t')
        if A > 0:
            d0 = d0[:A]+d0[(A+16):]
        if B > 0:
            d1 = d1[:B]+d1[(B+16):]
        if C > 0:
            v0 = v0[:C]+v0[(C+16):]
        if E > 0:
            v1 = v1[:E]+v1[(E+16):]
        if G > 0:
            a0 = a0[:G]+a0[(G+16):]
        if H > 0:
            a1 = a1[:H]+a1[(H+16):]

    display(Markdown("### Deslocamentos"))
    display(Markdown("$ \qquad x_p $ = $ %s $" %d0))
    display(Markdown("$ \qquad y_p $ = $ %s $" %d1))
    display(Markdown("### Velocidades"))
    display(Markdown("$ \qquad \dot{x_p} $ = $ %s $" %v0))
    display(Markdown("$ \qquad \dot{y_p} $ = $ %s $" %v1))
    display(Markdown("### Acelerações"))
    display(Markdown("$ \qquad \ddot{x_p} $ = $ %s $" %a0))
    display(Markdown("$ \qquad \ddot{y_p} $ = $ %s $" %a1))
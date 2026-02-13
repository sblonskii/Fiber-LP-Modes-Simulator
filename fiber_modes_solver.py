# Interaktywny symulator modów w światłowodzie włóknistym

import re
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from scipy.special import jv, kv
from scipy.optimize import brentq

mp.mp.dps = 50   #Ustawienie czułości obliczeń dla biblioteki mpmath -bardzo duża dokładność

"""
 Na początku ma miejsce parsowanie jednostek. Dzięki temu użytkownik otrzymuje dużą swobodę we wpisywaniu początkowych danych i może
to robić różnymi sposobami. Dodatkowo ułatwia to obliczenia dla komputera.
Objaśnienie: Parsowanie danych polega na tym, że komputer z ciągu tekstowego potrafi wyodrębnić np. liczby i jednostki.

"""
def Parsuj1(x):
    x = x.replace(" ","").replace("μ","u").lower()
    m = re.match(r"^\s*([+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?)\s*([a-z]*)\s*$", x)   #definicja przyjętego wzorca
    if not m:      #Jeśli podany ciąg tekstowy nie spełnia warunku to zwraca błąd
        raise ValueError(f"Nie rozpoznano wartości i jednostki wyrażenia '{x}'")
    wartosc = float(m.group(1))
    jednostka = m.group(2)
    if jednostka in ("","m"):
        return wartosc
    elif jednostka in ("um","mikrometr","micrometer","micron","u","mikrom"):
        return wartosc * 1e-6
    elif jednostka in ("mm"):
        return wartosc * 1e-3
    elif jednostka in ("nm"):
        return wartosc * 1e-9
    elif jednostka in ("pm"):
        return wartosc * 1e-12
    else:
        raise ValueError(f"Nieznana jednostka '{jednostka}'")

def Parsuj2(x):
    return float(x.strip())         # Pierwsza funkcja parsowania jest dla liczb w domyśle z jednostkami, a druga dla liczb bez jednostek

def ask_with_units(prompt: str, default_str: str, parser=Parsuj1):
    """Pobiera tekst z input(), akceptuje Enter=domyślne i parsuje jednostki."""
    raw = input(f"{prompt} [{default_str}]: ").strip()
    if raw == "":
        raw = default_str
    try:
        return parser(raw)
    except Exception as e:
        print(f"Błąd: {e}. Spróbuj ponownie.")
        return ask_with_units(prompt, default_str, parser)

def ask_int(prompt: str, default_val: int):
    raw = input(f"{prompt} [{default_val}]: ").strip()
    if raw == "":
        return default_val
    try:
        v = int(raw)
        return v
    except:
        print("Podaj liczbę całkowitą.")
        return ask_int(prompt, default_val)

# ====== INTERAKTYWNY INPUT ======
print("\n=== Parametry światłowodu ===\n")
print("Cześć! Jeśli chcesz skorzystać z programu, podaj proszę parametry,")
print("abym mógł wyznaczyć dla Ciebie rozkłady pól EM dla modów w światłowodzie.\n")
while True:
    try:
        lam = Parsuj1(input("Podaj długość fali λ (np. 1550 nm) : "))
        n1  = Parsuj2(input("Podaj współczynnik załamania rdzenia n1 (np. 1.450) : "))
        n2  = Parsuj2(input("Podaj współczynnik załamania płaszcza n2 (np. 1.442) : "))
        r   = Parsuj1(input("Podaj promień rdzenia r (np. 10 um) : "))

        # === WALIDACJA ===
        if n1 <= 1 or n2 <= 1:
            print("Współczynniki załamania muszą być większe od 1.\n")
            continue
        if n2 >= n1:
            print("Warunek prowadzenia nie spełniony: n1 musi być > n2.\n")
            continue
        break  # jeśli wszystko OK → wychodzimy z pętli

    except ValueError as e:
        print(f"Błąd: {e}\nSpróbuj ponownie.\n")
    except KeyboardInterrupt:
        print("\nProgram przerwany przez użytkownika (Ctrl+C).")
        exit()

# Po wyjściu z pętli wyswietla przyjęte dane:
print("\n Dane przyjęte:")
print(f"λ = {lam*1e9:.1f} nm")
print(f"n1 = {n1:.6f}, n2 = {n2:.6f}")
print(f"a = {r*1e6:.2f} µm")
print()


#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------



k0 = 2 * np.pi / lam
k1 = k0 * n1
k2 = k0 * n2

def U_od_beta(beta):
    return r * np.sqrt(k1**2 - beta**2)

def W_od_beta(beta):
    return r * np.sqrt(beta**2 - k2**2)

def fBessela_od_bety(l, beta):
    if not (k2< beta < k1):
        return np.nan
    u = U_od_beta(beta)
    w = W_od_beta(beta)
    if l == 0:
        rdzen = u * jv(-1, u) / jv(0, u)
        plaszcz = w * kv(-1, w) / kv(0, w)
    else:
        rdzen = u * jv(l-1, u) / jv(l, u)
        plaszcz = w * kv(l-1, w) / kv(l, w)

    return rdzen + plaszcz



Lmax = 15
for l in range(Lmax + 1):
    
    beta = np.linspace(k2*1.0001, k1*0.9999, 100000)
    
    wartosc = np.array([fBessela_od_bety(l, b) for b in beta])

    plt.plot(beta, wartosc, color = 'm', linewidth = 1)
    plt.axhline(y = 0, color = 'k', lw = 0.7)
    plt.xlabel("β")
    plt.ylabel("f(β)")
    plt.grid(True)
    plt.title("Wykres równania dyspersyjnego")
    plt.xlim(k2, k1)
    plt.ylim(-50, 50)




    miejsca_zerowe = []

    for i in range(len(beta)-1):
        y1 = wartosc[i]
        y2 = wartosc[i+1]
        if np.isfinite(y1) and np.isfinite(y2) and y1*y2 < 0 and y1 < 100 and y2 < 100  :
            try:
                x = brentq(lambda B: fBessela_od_bety(l,B), beta[i], beta[i+1])
                if not miejsca_zerowe or abs(x - miejsca_zerowe[-1]) > miejsca_zerowe[-1] * 1e-6:
                    miejsca_zerowe.append(x)
            except ValueError:
                pass
    if(len(miejsca_zerowe) > 0):
        if miejsca_zerowe:
            miejsca_zerowe = sorted(miejsca_zerowe, reverse=True)
        print("Znalezione β (mody LP dla l = {}) :".format(l))

        for p,b in enumerate(miejsca_zerowe, start = 1):
            print(f"p = {p} :  β = {b:.6e}   n_eff = {b/k0:.6f} ")



    def funkcja_r(beta, l, a=r, Rmax = 3.0, N = 1000):

        u = U_od_beta(beta)
        w = W_od_beta(beta)

        rr = np.linspace(0.0, Rmax*a, N)
        F = np.zeros_like(rr, dtype = float)

        mask_core = rr <= a
        rc = rr[mask_core]
        rl = rr[~mask_core]

        Fc = jv(l, u * rc / a)

        J_at_a = jv(l, u)
        K_at_a = kv(l, w) if w > 0 else 1.0
        B = J_at_a / K_at_a if np.isfinite(K_at_a) and K_at_a != 0 else 0.0

        Fl = B * kv(l, w * rl / a)

        F[mask_core] = Fc
        F[~mask_core] = Fl

        I = F**2
        m = np.nanmax(I)
        if np.isfinite(m) and m > 0:
            I = I / m

        return rr, I


    #Rysowanie profili radialnych dla znalezionych miejsc zerowych - bet dla których propagują się mody w światłowodzie o wskazanych parametrach.

    if len(miejsca_zerowe):

        for i, b, in enumerate(miejsca_zerowe, start = 1):
            if i<4:
                rr, F = funkcja_r(b, l, a = r, Rmax = 3.0)
                plt.figure()
                plt.plot(rr*1e6, F, color = 'k', lw = 1)
                plt.grid(True)
                plt.axhline(y = 0, color = 'k', ls = "--", lw = 0.7)
                plt.xlabel("r[μm]")
                plt.ylabel("Intensywność (unorm.)")
                plt.title(f"Profil radialny intensywności | l = {l} | p = {i} | n_eff ={b/k0:.5f}")


            else:
                pass

    def mapa_2d(beta, l, a = r, Rmax = 3.0, N = 1000):
        u = U_od_beta(beta)
        w = W_od_beta(beta)

        L = Rmax*a
        x = np.linspace(-L, L, N)
        y = np.linspace(-L, L, N)
        X,Y = np.meshgrid(x,y)
        R = np.hypot(X,Y)
        TH = np.arctan2(Y,X)


        F = np.zeros_like(R, dtype = float)
        mask_core = R <= a

        F[mask_core] = jv(l, u * R[mask_core] / a)
        B = jv(l, u) / kv(l, w) if w > 0 else 1.0
        if w > 0: F[~mask_core] = B * kv(l, w * R[~mask_core] / a)

       # zależność kątowa LP ~ cos(l*theta)
        ang = np.cos(l*TH) if l > 0 else 1.0
        I = (F * ang)**2
        return X, Y, I

    # przykład: pierwsza znaleziona beta
    if len(miejsca_zerowe):
        for i, b in enumerate(miejsca_zerowe, start=1):
            if i < 20:
                n_eff = b/k0
                if n_eff > n2:
                    X, Y, I = mapa_2d(b, l, a=r)
                    plt.figure()
                    plt.pcolormesh(X*1e6, Y*1e6, I, shading='auto')
                    c = plt.Circle((0,0), r*1e6, color='w', fill=False, ls='--', lw=0.8)
                    plt.gca().add_artist(c); plt.gca().set_aspect('equal')
                    plt.xlabel("x [µm]")
                    plt.ylabel("y [µm]")
                    plt.title(f"Intensywność 2D | l={l}, p={i}")
                    plt.tight_layout()
                    plt.show()
                else:
                    pass
            else:
                pass

print()
print()
print("Dzięki za skorzystanie z programu!")

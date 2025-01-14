
# Modèles mathématiques 1D et 2D pour l’écoulement sanguin

## Modèle 1D
Le modèle 1D repose sur une **intégration par section des équations de Navier-Stokes** sous l’hypothèse d’un **écoulement incompressible** dans des artères supposées fines. Ce modèle est particulièrement adapté pour des études globales du réseau artériel, où la géométrie est approximativement linéaire ou faiblement courbée.

### Hypothèses et simplifications
- L’écoulement est considéré comme **incompressible**.
- L’artère est modélisée comme un tube cylindrique de section variable en fonction de la pression.
- Un profil de vitesse parabolique est utilisé, permettant une **moyennisation** sur la section transversale de l’artère.

### Équations principales
Les équations dérivées sont un système d’**équations hyperboliques** aux dérivées partielles décrivant la conservation de la masse et de la quantité de mouvement :

1. **Conservation de la masse** :
```math
∂_t A + ∂_x Q = 0
```

2. **Conservation de la quantité de mouvement** :
```math
∂_t Q + ∂_x \left( \alpha \frac{Q^2}{A} + \frac{1}{\rho} A P(A, x) \right) - \partial_x \left( 3\nu A \partial_x\left(\frac Q A\right) \right) = \frac{1}{\rho} P(A, x) ∂_x A - \frac{2\pi R K}{1-\frac{Rk}{4\nu}} \frac{Q}{A}
```

### Énergie et relation d’entropie du modèle 1D
L’énergie associée au système est donnée par :
```math
E(t, x) = \frac{A u_x^2}{2} + \frac{1}{\rho} A P(A, x) - \frac{\beta(x)}{3 \rho A_0(x)} A^{3/2}
```

La relation d’entropie vérifiée par cette énergie est :
```math
∂_t E + ∂_x \left( \left( E + \frac{\beta(x)}{3 \rho A_0(x)} A^{3/2} \right) u_x \right) = ∂_x \left( 3 \nu A ∂_x \left( \frac{Q}{A} \right) \right) u_x + \frac{2 \pi R k}{1 - R k / 4 \nu} u_x^2 ≤ 0
```

Sous des conditions aux limites nulles :
```math
∂_t \left( \int_0^L E \, dx \right) = - 3 \nu \int_0^L A (∂_x u_x)^2 \, dx - \frac{2 \pi R k}{1 - R k / 4 \nu} \int_0^L u_x^2 \, dx < 0
```

---

## Modèle 2D
Le modèle 2D est dérivé à partir d’une **intégration radiale des équations de Navier-Stokes**, permettant de mieux représenter les effets locaux dans des configurations géométriques complexes, comme les **bifurcations artérielles** et les **anévrismes sévères**.

### Hypothèses et simplifications
- L’écoulement est supposé **incompressible**.
- La géométrie de l’artère est décrite à l’aide d’un système de coordonnées curvilignes (\( s, 	heta \)).
- Le profil de vitesse est obtenu sans recourir à un ansatz spécifique.

### Équations principales
1. **Conservation de la masse** :
```math
∂_t A + ∂_θ \left( \frac{Q_{Rθ}}{A} \right) + ∂_s(Q_s) = 0
```

2. **Conservation de la quantité de mouvement (composante radiale et axiale)** :
```math
∂_t (Q_{Rθ}) + ∂_θ \left( \frac{Q_{Rθ}^2}{2 A^2} + A P \right) + ∂_s \left( \frac{Q_{Rθ} Q_s}{A} \right) = \frac{2 R}{3} C \sin θ \frac{Q_s^2}{A} + \frac{2 R k Q_{Rθ}}{A} + P∂_θ (A)
```
```math
∂_t (Q_s) + ∂_θ \left( \frac{Q_s Q_{Rθ}}{A^2} \right) + ∂_s \left( \frac{Q_s^2}{A} - \frac{Q_{Rθ}^2}{2 A^2} + A P \right) = - \frac{2 R}{3} C \sin θ \frac{Q_{Rθ} Q_s}{A^2} + \frac{k R Q_s}{A} + P∂_s (A)
```

### Énergie et relation d’entropie du modèle 2D
L’énergie associée au système est donnée par :
```math
E(t, θ, s) = A \left( \frac{9}{8} u_θ^2 + \frac{u_s^2}{2} + p \right) - p̃
```

La relation d’entropie correspondante est :
```math
∂_t E + ∂_θ \left( \frac{3}{2} \frac{u_θ}{R} \left( E + p̃ - \frac{9}{16} A u_θ^2 \right) \right) + ∂_s \left( u_s \left( E + p̃ - \frac{9}{16} A u_θ^2 \right) \right) = \frac{9}{4} R k u_θ^2 + k R u_s^2 ≤ 0
```

Cette relation garantit que l’énergie décroît localement dans le temps, ce qui assure la **stabilité** du modèle.

---

## Comparaison des modèles 1D et 2D
- **Modèle 1D** :
  - Rapide et efficace pour des simulations globales sur de grands réseaux artériels.
  - Bien adapté pour des géométries simples ou faiblement courbées.
  - Coût de calcul très faible.
- **Modèle 2D** :
  - Plus précis pour des géométries complexes (bifurcations, anévrismes).
  - Permet de mieux capturer les effets locaux et les interactions fluide-structure.
  - Coût de calcul modéré par rapport aux modèles tridimensionnels (NS-FSI 3D).

L’utilisation combinée de ces deux modèles permet une **alternative efficace aux simulations 3D**, tout en offrant un bon compromis entre précision et coût de calcul.

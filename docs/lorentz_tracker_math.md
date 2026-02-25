# Lorentz Tracker Algebra (Heuristic And SU2-Phase-Aware View)

This note defines the exact matrix conventions used by the tracker and summarizes how relative Wigner angles are obtained.

## 1. Basis and Action Convention

We use column four-vectors in the basis
\[
v = (p_x, p_y, p_z, E)^\top.
\]

A tracked Lorentz map acts as
\[
v' = \Lambda v.
\]

If an instruction path has steps \(\Lambda_1,\Lambda_2,\ldots,\Lambda_n\), the accumulated map is left-composed:
\[
\Lambda_{\mathrm{tot}} = \Lambda_n \cdots \Lambda_2 \Lambda_1.
\]

For two paths (reference \(r\), other \(o\)), the relative map is
\[
\Delta = \Lambda_o \Lambda_r^{-1}.
\]

This is exactly `relative = other * inv(reference)`.

## 2. Elementary 4x4 Matrices Used

In \((p_x,p_y,p_z,E)\), the active rotations/boost are
\[
R_z(\phi)=
\begin{pmatrix}
\cos\phi & -\sin\phi & 0 & 0 \\
\sin\phi & \cos\phi  & 0 & 0 \\
0&0&1&0\\
0&0&0&1
\end{pmatrix},
\]
\[
R_y(\theta)=
\begin{pmatrix}
\cos\theta & 0 & \sin\theta & 0 \\
0&1&0&0\\
-\sin\theta & 0 & \cos\theta & 0 \\
0&0&0&1
\end{pmatrix},
\]
\[
B_z(\xi)=
\begin{pmatrix}
1&0&0&0\\
0&1&0&0\\
0&0&\cosh\xi&\sinh\xi\\
0&0&\sinh\xi&\cosh\xi
\end{pmatrix}.
\]

## 3. Lorentz Decode to Helicity Parameters

Given \(M\), boost parameters are decoded from the 4th column \(M[:,4]\) (image of rest vector \((0,0,0,1)^\top\)):
\[
\gamma = M_{44},\qquad
\xi = \operatorname{arcosh}(\gamma),
\]
\[
\phi = \operatorname{atan2}(M_{24},M_{14}),\qquad
\theta = \arccos\!\left(\frac{M_{34}}{\sqrt{M_{14}^2+M_{24}^2+M_{34}^2}}\right).
\]

Then remove the boost part:
\[
M_{\mathrm{rf}} = B_z(-\xi)\,R_y(-\theta)\,R_z(-\phi)\,M.
\]

From the \(3\times3\) rotation block \(R\) of \(M_{\mathrm{rf}}\), decode ZYZ:
\[
\phi_{\mathrm{rf}}=\operatorname{atan2}(R_{23},R_{13}),\quad
\theta_{\mathrm{rf}}=\arccos(R_{33}),\quad
\psi_{\mathrm{rf}}=\operatorname{atan2}(R_{32},-R_{31}).
\]

For phase-aware convention, we normalize
\[
\psi_{\mathrm{rf}}\in[-\pi,3\pi)
\]
via
\[
\psi \mapsto \operatorname{mod}(\psi+\pi,4\pi)-\pi.
\]

## 4. SU2 Companion Tracking

The tracker also accumulates a \(2\times2\) SU(2)-like matrix \(U\) (spinor representation):
\[
U_{\mathrm{tot}} = U_n \cdots U_2 U_1.
\]

Elementary factors:
\[
U_{R_z}(\phi)=
\begin{pmatrix}
e^{-i\phi/2}&0\\
0&e^{i\phi/2}
\end{pmatrix},\quad
U_{R_y}(\theta)=
\begin{pmatrix}
\cos(\theta/2)&-\sin(\theta/2)\\
\sin(\theta/2)&\cos(\theta/2)
\end{pmatrix},
\]
\[
U_{B_z}(\xi)=
\begin{pmatrix}
e^{\xi/2}&0\\
0&e^{-\xi/2}
\end{pmatrix}.
\]

The key point is covering:
\[
\mathrm{SU}(2)\to \mathrm{SO}(3)
\]
is 2-to-1, so \(U\) and \(-U\) map to the same spatial rotation. Therefore SO(3)/4x4-only decoding can miss a \(2\pi\) spinor phase branch, while SU2 tracking can preserve it.

## 5. How SU2 Resolves the Missing \(2\pi\) Branch

Write Euler part as
\[
U_{\mathrm{rot}}=U_{R_z}(\phi_{\mathrm{rf}})\,U_{R_y}(\theta_{\mathrm{rf}})\,U_{R_z}(\psi_{\mathrm{rf}}).
\]

Because \(U\) is single-valued over \(4\pi\), \(\psi_{\mathrm{rf}}\) can be selected on a \(4\pi\)-length interval (here \([-\pi,3\pi)\)) consistently, instead of modulo \(2\pi\) only.

In the current implementation:
- SU2 is accumulated in the tracker,
- for pure rotations (\(\xi \approx 0\)), \(U\) is used to choose the \(\psi\) branch (\(\psi\) vs \(\psi+2\pi\)),
- for generic boosted transforms, decoding remains \(\Lambda\)-branch based.

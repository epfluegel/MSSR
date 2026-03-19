# MSSR Matrix Share Size Reduction Implementations

This repository contains Maple reference implementations of the algorithms studied in the paper *Online Secret Sharing Schemes Based on Matrix Normal Forms*. The idea is to represent a secret as a matrix over a finite field and then use similarity transformations to separate it into a public component and a smaller private component. In the probabilistic MSSR construction, the algorithm searches for a cyclic-vector representation, with expected fast termination (less than three iterations in the paper’s analysis) and reduced private share size on the order of \(\Theta(\sqrt{|s|})\). This gives practical short-share behaviour while preserving the intended reconstruction guarantees.

The repository also includes a deterministic variant based on the Frobenius canonical normal form, which avoids probabilistic cyclic-vector selection and applies to every input matrix. This deterministic version has cubic matrix-algebra complexity (\(O(t^3)\), with \(t\) the matrix dimension), making it a predictable baseline when guaranteed canonicalisation is required. For first-time users, a useful way to read the code is to follow the same approach as in the paper: secret-to-matrix conversion over \(\mathbb{F}_p\), normal-form/similarity computation, extraction of compact invariants, and finally share generation. This highlights the trade-off between speed (probabilistic method), guaranteed applicability (deterministic method), and security instantiation (e.g., ramp-style sharing).

## Getting Started

If this is your first time in the repository, use the workflow below to connect the implementation to the paper's theory.

1. **Choose the algebraic parameters**: start from a secret \(s\), select a prime \(p\), and map \(s\) into a matrix over \(\mathbb{F}_p\) with dimension \(t\).
2. **Pick an algorithmic path**: use the probabilistic MSSR path when you want fast expected runtime, or the deterministic Frobenius path when you need guaranteed canonicalisation.
3. **Compute the normal-form representation**: obtain either a companion/cyclic form or a Frobenius block form via similarity transformation.
4. **Extract compact private data**: read the invariant coefficients from the transformed matrix to produce the shortened private payload.
5. **Generate shares**: feed the private payload into your secret-sharing backend (e.g., threshold or ramp configuration), while storing the public component as authenticated public data.

A practical way to validate understanding is to run both paths on the same \((s,p)\) input and compare runtime, output size, and reconstruction correctness.

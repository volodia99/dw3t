import os

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

plt.style.use("nonos.default")

def read_radmc3d_opacity(*, filename:str, directory:str|None=None):
    """
        Adapted from fargo2radmc3d
    """
    if directory is None:
        directory = "."
    path = os.path.join(directory, filename)
    head = 3
    header = np.fromfile(path, dtype=int, count=head, sep=" ")
    Nlam = header[1]
    lbda, kappa_abs, kappa_sca, g = np.fromfile(
        path,
        dtype=float,
        sep=" ",
    )[head:head+4*Nlam].reshape(Nlam, 4).T
    return (lbda, kappa_abs, kappa_sca, g)

def plot_radmc3d_opacity(*, dustsize:np.ndarray[tuple[int], np.dtype[u.Quantity]], at_lambda:u.Quantity, directory:str|None=None):
    """
        Adapted from fargo2radmc3d
    """
    if directory is None:
        directory = "."
    fig, ax = plt.subplots(figsize=(8,5))
    ax.set(
        xscale="log",
        yscale="log",
        xlabel=f"Dust size [{dustsize.unit}]",
        ylabel=r"Opacities $[{\rm cm}^2.{\rm g}^{-1}]$",
    )

    Nspec = dustsize.shape[0]
    # computes the 'magnitude' of Nspec
    mag = int(np.ceil(np.log10(Nspec)))

    absorption = np.zeros(Nspec)
    scattering = np.zeros(Nspec)
    for kk in range(Nspec):
        filename = f"dustkapscatmat_{str(kk).zfill(mag)}.inp"
        lbda, kappa_abs, kappa_sca, g = read_radmc3d_opacity(filename=filename, directory=directory)

        lbda1 = at_lambda.to(u.micron).value
        i1 = np.argmin(np.abs(lbda-lbda1))

        # interpolation in log
        #TODO: to be improved
        l1 = lbda[i1-1]
        l2 = lbda[i1+1]
        k1 = kappa_abs[i1-1]
        k2 = kappa_abs[i1+1]
        ks1 = kappa_sca[i1-1]
        ks2 = kappa_sca[i1+1]
        absorption[kk] =  (k1*np.log(l2/lbda1) +  k2*np.log(lbda1/l1))/np.log(l2/l1)
        scattering[kk] = (ks1*np.log(l2/lbda1) + ks2*np.log(lbda1/l1))/np.log(l2/l1)

    ax.plot(
        dustsize, 
        absorption, 
        lw=2., 
        linestyle="-", 
        color="k", 
        label=rf"$\kappa_{{abs}}$ at {at_lambda.to(u.mm)}",
    )
    ax.plot(
        dustsize, 
        absorption+scattering, 
        lw=2., 
        linestyle="--", 
        color="k", 
        label=rf"$\kappa_{{abs}}+\kappa_{{sca}}$ at {at_lambda.to(u.mm)}",
    )

    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(directory, "opacities.pdf"), dpi=150)

# import numpy as np
# import matplotlib.pyplot as plt
# plot_opacities(amin=1e-8,amax=1e-3,nbin=30,lbda1=870)

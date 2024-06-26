(sec-install-windows)=
# Windows Packages

Windows installers are provided for stable versions of Cantera. These installers
provide header/library files that can be used to compile C++ applications.

:::{attention}
The *legacy* Matlab Cantera interface is discontinued and removed in Cantera 3.1. Users
requiring support of legacy Matlab Cantera code should continue using Cantera 3.0
packages, or migrate their code base to the experimental Matlab toolbox that is
currently under development.
:::

:::{seealso}
To install the Cantera Python package, see the [pip](pip) or [conda](conda)
installation instructions. The Python package is required if:

- You need to convert Chemkin-format input files to YAML
- You need to convert legacy CTI or XML input files to YAML
:::

1. **Remove old versions of Cantera**

- Use The Windows "Add/Remove Programs" interface to remove previous versions of
  the `Cantera` package.

2. **Install Cantera**

- Go to the [Cantera Releases](https://github.com/Cantera/cantera/releases) page and
  download **Cantera-3.0.0-x64.msi**.
- Run the installer and follow the prompts.

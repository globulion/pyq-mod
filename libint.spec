Name:		libint
Version:	1.1.4
Release:	4%{?dist}
Summary:	A library for computing electron repulsion integrals efficiently
Group:		System Environment/Libraries
License:	GPLv2+
URL:		http://www.files.chem.vt.edu/chem-dept/valeev/software/libint/libint.html
Source0:	http://www.files.chem.vt.edu/chem-dept/valeev/software/libint/src/libint-%{version}.tar.gz
BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)
# Required to build documentation
BuildRequires:	/usr/bin/bibtex
BuildRequires:	/usr/bin/pdflatex

# g++ crashes on pcc when trying to build in EPEL-5
ExcludeArch:	ppc

%description
LIBINT computes the Coulomb and exchange integrals, which in electronic
structure theory are called electron repulsion integrals (ERIs). This is by
far the most common type of integrals in molecular structure theory.

LIBINT uses recursive schemes that originate in seminal Obara-Saika method and
Head-Gordon and Pople’s variation thereof. The idea of LIBINT is to optimize
computer implementation of such methods by implementing an optimizing compiler
to generate automatically highly-specialized code that runs well on
super-scalar architectures.

%package devel
Summary:	Development headers and libraries for libint
Group:		Development/Libraries
Requires:	libint = %{version}-%{release}
Requires:	libderiv = %{version}-%{release}
Requires:	libr12 = %{version}-%{release}

%description devel
This package contains development headers and libraries for libint.
It also contains a programmer's manual.

%package -n libr12
Summary:	A library for computing integrals that arise in Kutzelnigg’s linear R12 theories
Group:		System Environment/Libraries

%description -n libr12
libr12 computes types integrals that appear in Kutzelnigg’s linear R12 theories
for electronic structure. All linear R12 methods, such as MP2-R12, contain
terms in the wave function that are linear in the inter-electronic distances
r_{ij} (hence the name). Appearance of several types of two-body integrals is
due to the use of the approximate resolution of the identity to reduce three-
and four-body integrals to products of simpler integrals.

%package -n libderiv
Summary:	A library for computing derivatives of electron repulsion integrals
Group:		System Environment/Libraries

%description -n libderiv
libderiv computes first and second derivatives of ERIs with respect to the
coordinates of the basis function origin. This type of integrals are also very
common in electronic structure theory, where they appear in analytic gradient
expressions. The derivatives are typically used in the calculation of forces.


%prep
%setup -q

%build
%configure --enable-shared --disable-static \
 --with-libint-max-am=6 --with-libderiv-max-am1=5 --with-libderiv-max-am2=4 \
 --with-libr12-max-am=5 

%ifarch s390
# change to -O1 to prevent memory exhaustion by g++
%global optflags %(echo %{optflags} | sed 's/-O2/-O1/')
%endif

make CFLAGS="%{optflags}" CXXFLAGS="%{optflags}" %{?_smp_mflags}

# Build documentation
cd doc/progman
pdflatex progman
bibtex progman
pdflatex progman
pdflatex progman


%install
rm -rf %{buildroot} 
make install DESTDIR=%{buildroot}
find %{buildroot} -name *.la -delete
find %{buildroot} -name *.so.*.* -exec chmod 755 {} \;

%clean
rm -rf %{buildroot}

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%post -n libderiv -p /sbin/ldconfig
%postun -n libderiv -p /sbin/ldconfig

%post -n libr12 -p /sbin/ldconfig
%postun -n libr12 -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%doc LICENSE
%{_libdir}/libint*.so.*

%files -n libderiv
%defattr(-,root,root,-)
%{_libdir}/libderiv*.so.*

%files -n libr12
%defattr(-,root,root,-)
%{_libdir}/libr12*.so.*

%files devel
%defattr(-,root,root,-)
%doc doc/progman/progman.pdf
%{_includedir}/libint/
%{_includedir}/libderiv/
%{_includedir}/libr12/
%{_libdir}/*.so


%changelog
* Mon Dec 13 2010 Dan Horák <dan[at]danny.cz> - 1.1.4-4
- workaround memory exhaustion on s390

* Tue Nov 30 2010 Jussi Lehtola <jussilehtola@fedoraproject.org> - 1.1.4-3
- Increase maximum angular momentum values by 2, making it possible to
  use basis sets that use up to I-type functions, such as Dunning's cc-pVXZ
  basis sets.
- Split libderiv and libr12 into their own packages, as e.g. PyQuante currently
  only needs the libint library.

* Fri Jul 24 2009 Fedora Release Engineering <rel-eng@lists.fedoraproject.org> - 1.1.4-2
- Rebuilt for https://fedoraproject.org/wiki/Fedora_12_Mass_Rebuild

* Tue May 26 2009 Jussi Lehtola <jussilehtola@fedoraproject.org> - 1.1.4-1
- First release.

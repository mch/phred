# phred.spec
#
# Copyright (c) 2004 Matt Hughes mhughe@uvic.ca
#
%define name phred
%define version 0.1.2
%define release 1
%define manifest %{_builddir}/%{name}-%{version}-%{release}.manifest

# required items
Name: %{name}
Version: %{version}
Release: %{release}
Copyright: GPL
Group: Application/Misc

# optional items
#Vendor: Matt Hughes
#Distribution:
#Icon:
#URL:
#Packager: Matt Hughes mhughe@uvic.ca

# source + patches
Source: %{name}-%{version}.tar.gz
#Source1:
#Patch:
#Patch1:

# RPM info
#Provides:
#Requires:
#Conflicts:
#Prereq:

#Prefix: /usr
BuildRoot: /var/tmp/%{name}-%{version}

Summary: Phred is a parallel finite difference time domain electromagnetics simulator.

%description
Phred is a parallel finite difference time domain electromagnetics simulator.
phred.spec.  Generated from phred.spec.in by configure.
Please edit phred.spec.in to add several more lines of description
here if appropriate, and to delete these instructions.

%prep
%setup -q
#%patch0 -p1

%build
%configure
make

%install
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT
%makeinstall

# __os_install_post is implicitly expanded after the
# %install section... do it now, and then disable it,
# so all work is done before building manifest.

%{?__os_install_post}
%define __os_install_post %{nil}

# build the file list automagically into %{manifest}

cd $RPM_BUILD_ROOT
rm -f %{manifest}
find . -type d \
        | sed '1,2d;s,^\.,\%attr(-\,root\,root) \%dir ,' >> %{manifest}
find . -type f \
        | sed 's,^\.,\%attr(-\,root\,root) ,' >> %{manifest}
find . -type l \
        | sed 's,^\.,\%attr(-\,root\,root) ,' >> %{manifest}

#%pre
#%post
#%preun
#%postun

%clean
rm -f %{manifest}
rm -rf $RPM_BUILD_ROOT

%files -f %{manifest}
%defattr(-,root,root)
#%doc README
#%docdir
#%config

%changelog

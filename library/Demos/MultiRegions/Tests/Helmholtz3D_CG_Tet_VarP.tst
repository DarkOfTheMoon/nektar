<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG for Tet with variable P</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Tet_VarP.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Tet_VarP.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.000703024</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.00601347</value>
        </metric>
    </metrics>
</test>


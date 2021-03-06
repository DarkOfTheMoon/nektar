<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interp field to a plane of points (also calculate cp and cp0)</description>
    <executable>FieldConvert</executable>
    <parameters>-e  -m interppoints:cp=0,0.5:plane=10,10,0.1,-0.9,-0.9,0.1,0.9,-0.9,0.1,0.9,0.9,0.1,-0.9,0.9:fromxml=chan3D.xml:fromfld=chan3D.fld chan3D_plane.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6"> 0.1</value>
            <value variable="y" tolerance="1e-6">0.574456</value>
            <value variable="z" tolerance="1e-6">0.574456</value>
            <value variable="u" tolerance="1e-6">0.730329</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">1.8</value>
            <value variable="Cp" tolerance="1e-6">3.6</value>
            <value variable="Cp0" tolerance="1e-6">4.14809</value>
        </metric>
    </metrics>
</test>


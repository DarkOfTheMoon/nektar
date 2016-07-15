<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_SubStep.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_SubStep.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.000450925</value>
            <value variable="v" tolerance="1e-6">0.00015642</value>
            <value variable="w" tolerance="1e-6">7.01706e-05</value>
            <value variable="p" tolerance="1e-6">0.000609814</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00929358</value>
            <value variable="v" tolerance="1e-6">0.0014178</value>
            <value variable="w" tolerance="1e-6">0.000511184</value>
            <value variable="p" tolerance="2e-6">0.00855209</value>
        </metric>
    </metrics>
</test>

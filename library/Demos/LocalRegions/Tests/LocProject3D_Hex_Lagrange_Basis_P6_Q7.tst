<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Hex Lagrange basis P=6 Q=7</description>
    <executable>LocProject3D</executable>
    <parameters>8 10 10 10 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">2.04616e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">8.52651e-13</value>
        </metric>
    </metrics>
</test>

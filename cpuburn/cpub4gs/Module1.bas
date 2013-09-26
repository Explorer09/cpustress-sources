Attribute VB_Name = "Module1"

Public Declare Function StartBurn Lib "CPUburn.DLL" (ByRef whichone As Long, ByRef priority As Long, mem_siz As Long) As Long
Public Declare Function Kill_me Lib "CPUburn.DLL" () As Long

Public status As Long, tp As Long, pr As Long, ms As Long



VERSION 5.00
Object = "{831FDD16-0C5C-11D2-A9FC-0000F8754DA1}#2.0#0"; "MSCOMCTL.OCX"
Begin VB.Form Form1 
   BorderStyle     =   4  'Fixed ToolWindow
   Caption         =   "CPUburn"
   ClientHeight    =   3156
   ClientLeft      =   36
   ClientTop       =   276
   ClientWidth     =   3960
   Icon            =   "Form1.frx":0000
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   NegotiateMenus  =   0   'False
   ScaleHeight     =   3156
   ScaleWidth      =   3960
   StartUpPosition =   3  'Windows Default
   Begin VB.OptionButton Option1 
      Caption         =   "K7"
      Height          =   252
      Index           =   6
      Left            =   360
      TabIndex        =   15
      Top             =   1800
      Width           =   612
   End
   Begin VB.OptionButton Option1 
      Caption         =   "BX"
      Height          =   252
      Index           =   5
      Left            =   360
      TabIndex        =   14
      Top             =   1560
      Width           =   612
   End
   Begin VB.OptionButton Option1 
      Caption         =   "MMX"
      Height          =   252
      Index           =   4
      Left            =   360
      TabIndex        =   13
      Top             =   1320
      Width           =   800
   End
   Begin VB.OptionButton Option1 
      Caption         =   "P5"
      Height          =   252
      Index           =   3
      Left            =   360
      TabIndex        =   12
      Top             =   1080
      Width           =   612
   End
   Begin VB.OptionButton Option1 
      Caption         =   "K6"
      Height          =   252
      Index           =   2
      Left            =   360
      TabIndex        =   11
      Top             =   840
      Width           =   612
   End
   Begin VB.OptionButton Option1 
      Caption         =   "P6"
      Height          =   252
      Index           =   1
      Left            =   360
      TabIndex        =   10
      Top             =   600
      Value           =   -1  'True
      Width           =   612
   End
   Begin VB.CommandButton Command1 
      Caption         =   "BURN !!!"
      Height          =   495
      Left            =   2160
      TabIndex        =   6
      Top             =   2520
      Width           =   1695
   End
   Begin VB.Frame Frame2 
      Caption         =   "Process Priority:"
      Height          =   2292
      Left            =   2160
      TabIndex        =   1
      Top             =   120
      Width           =   1695
      Begin VB.OptionButton Option9 
         Caption         =   "Realtime"
         Enabled         =   0   'False
         Height          =   252
         Left            =   240
         TabIndex        =   5
         Top             =   1800
         Width           =   1190
      End
      Begin VB.OptionButton Option8 
         Caption         =   "High"
         Height          =   252
         Left            =   240
         TabIndex        =   4
         Top             =   1320
         Width           =   1095
      End
      Begin VB.OptionButton Option7 
         Caption         =   "Normal"
         Height          =   252
         Left            =   240
         TabIndex        =   3
         Top             =   840
         Width           =   1095
      End
      Begin VB.OptionButton Option6 
         Caption         =   "Idle"
         Height          =   252
         Left            =   240
         TabIndex        =   2
         Top             =   360
         Value           =   -1  'True
         Width           =   852
      End
   End
   Begin VB.Frame Frame1 
      Caption         =   "Test Type:"
      Height          =   2895
      Left            =   120
      TabIndex        =   0
      Top             =   120
      Width           =   1932
      Begin MSComctlLib.Slider Slider1 
         Height          =   1692
         Left            =   1080
         TabIndex        =   7
         ToolTipText     =   "Set Test RAM Size"
         Top             =   360
         Width           =   504
         _ExtentX        =   910
         _ExtentY        =   2985
         _Version        =   393216
         BorderStyle     =   1
         Enabled         =   0   'False
         Orientation     =   1
         Max             =   15
         TextPosition    =   1
      End
      Begin VB.Label Label2 
         Caption         =   " 64 kB"
         Height          =   252
         Left            =   1200
         TabIndex        =   9
         Top             =   2460
         Width           =   612
      End
      Begin VB.Label Label1 
         Caption         =   "BX && MMX RAM Size:"
         Height          =   600
         Left            =   120
         TabIndex        =   8
         Top             =   2220
         Width           =   1095
      End
   End
End
Attribute VB_Name = "Form1"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False


Private Sub Command1_Click()
Dim result As Long

If Option1(1).Value = True Then tp = 1
If Option1(2).Value = True Then tp = 2
If Option1(3).Value = True Then tp = 3
If Option1(4).Value = True Then tp = 4
If Option1(5).Value = True Then tp = 5
If Option1(6).Value = True Then tp = 6

If Option6.Value = True Then pr = 2
If Option7.Value = True Then pr = 1
If Option8.Value = True Then pr = 4
If Option9.Value = True Then pr = 5

If status = 0 Then
result = StartBurn(tp, pr, ms)
Command1.Caption = "STOP"
Else:
result = Kill_me()
Command1.Caption = "BURN !!!"
End If

status = status + 1
status = status And 1
End Sub

Private Sub Command2_Click()
frmAbout.Show

End Sub

Private Sub Form_Load()

If Command() = "-r" Then Option9.Enabled = True

status = 0
Slider1.Value = 5
ms = 5
End Sub

Private Sub Form_Terminate()
Dim r As Long

r = Kill_me()

End Sub





Private Sub Option1_Click(Index As Integer)
Slider1.Enabled = False
If Index = 4 Then Slider1.Enabled = True
If Index = 5 Then Slider1.Enabled = True
End Sub

Private Sub Option9_Click()
Dim dummy As String
dummy = MsgBox("Enabling this option will most likely" + Chr(13) + "make your computer unresponsive !", vbExclamation, "Warning !", "", 0)

End Sub

Private Sub StatusBar1_PanelClick(ByVal Panel As MSComctlLib.Panel)

End Sub


Private Sub Slider1_Scroll()
Dim st As String
Dim mem_size As Long

ms = Slider1.Value
mem_size = 2 ^ (Slider1.Value + 1)
If (mem_size) >= 1024 Then
mem_size = mem_size / 1024
st = Str(mem_size) + " MB"
Else
st = Str(mem_size) + " kB"
End If

Slider1.Text = st
Label2.Caption = st

End Sub

// CPUburn.cpp : Defines the entry point for the DLL application.
//

#include "stdafx.h"

int _stdcall StartBurn (DWORD,DWORD,DWORD);
int _stdcall IsStillBurning ();
int _stdcall Kill_me ();

char *burn_exe_name[100];
LPTSTR  command_line_string="\0     ";
STARTUPINFO startup_info;
PROCESS_INFORMATION proc_info;
DWORD dir_name_len;
char *current_dir[100];
char *const_names="\\BURNP6 \\BURNK6 \\BURNP5 \\BURNMX \\BURNBX \\BURNK7";
char *exe=".exe\0";
BYTE flags;

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}


int _stdcall StartBurn (DWORD whichone,DWORD priority,DWORD mem_size)
{

command_line_string="\0";
dir_name_len=GetCurrentDirectory(100,(char*)current_dir);

_asm{
	lea		edi,current_dir
	mov		edx,dir_name_len
	add		edi,edx
	mov		ebx,whichone
	mov		eax,dword ptr[ebx]
	push	eax
	dec		eax
	shl		eax,3
	mov		esi,const_names
	add		esi,eax
	mov		ecx,7
	rep	movsb
	add		ecx,5
	mov		esi,exe
	rep movsb
	mov		ecx,edx
	add		ecx,11
	lea		esi,current_dir
	lea		edi,burn_exe_name
	rep movsb

		mov		byte ptr [edi],32
		inc		edi
		pop		ebx
		cmp		ebx,4
		jl		skip_cls
		cmp		ebx,5
		jg		skip_cls
		//mov		edi,command_line_string
		mov		edx,mem_size	
		mov		eax,[edx]
		add		eax,65
		mov		byte ptr[edi],al
		mov		byte ptr[edi+1],0
skip_cls:	
}


startup_info.cb=sizeof(STARTUPINFO);
startup_info.lpReserved=NULL;
startup_info.lpDesktop="";
startup_info.lpTitle="Burning!";
startup_info.lpReserved2=NULL;
startup_info.cbReserved2=0;
startup_info.dwFlags=STARTF_USESHOWWINDOW;
startup_info.wShowWindow=SW_MINIMIZE;

_asm{
	mov		ebx,priority
	mov		cl,byte ptr[ebx]
	mov		edx,0x20
	dec		ecx
	shl		edx,cl
	mov		eax,edx
	and		eax,0x100
	shr		eax,1
	or		edx,eax
	mov		flags,dl


}
/*
HIGH_PRIORITY_CLASS = 0x80
IDLE_PRIORITY_CLASS = 0x40
NORMAL_PRIORITY_CLASS = 0x20
REALTIME_PRIORITY_CLASS = 0x00
*/

CreateProcess(NULL,burn_exe_name,NULL,NULL,TRUE,flags,
			  NULL,NULL,&startup_info,&proc_info);


return(0);
}


int _stdcall IsStillBurning ()
{


return(0);
}


int _stdcall Kill_me ()
{

TerminateProcess(proc_info.hProcess,0);

return(0);
}
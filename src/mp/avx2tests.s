	.file	"avx2tests.cpp"
	.text
	.type	_ZL5loadaPdi, @function
_ZL5loadaPdi:
.LFB850:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	movq	%rdi, -40(%rsp)
	movl	%esi, -44(%rsp)
	movl	-44(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rsp), %rax
	addq	%rdx, %rax
	movq	%rax, -24(%rsp)
	movq	-24(%rsp), %rax
	vmovapd	(%rax), %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE850:
	.size	_ZL5loadaPdi, .-_ZL5loadaPdi
	.type	_ZL5loaduPdi, @function
_ZL5loaduPdi:
.LFB852:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	movq	%rdi, -40(%rsp)
	movl	%esi, -44(%rsp)
	movl	-44(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rsp), %rax
	addq	%rdx, %rax
	movq	%rax, -24(%rsp)
	movq	-24(%rsp), %rax
	vmovupd	(%rax), %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE852:
	.size	_ZL5loaduPdi, .-_ZL5loaduPdi
	.type	_ZL6storeaU8__vectordPdi, @function
_ZL6storeaU8__vectordPdi:
.LFB854:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$12, %rsp
	vmovapd	%ymm0, -108(%rsp)
	movq	%rdi, -116(%rsp)
	movl	%esi, -120(%rsp)
	movl	-120(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-116(%rsp), %rax
	addq	%rdx, %rax
	movq	%rax, -20(%rsp)
	vmovapd	-108(%rsp), %ymm0
	vmovapd	%ymm0, -76(%rsp)
	movq	-20(%rsp), %rax
	vmovapd	-76(%rsp), %ymm0
	vmovapd	%ymm0, (%rax)
	nop
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE854:
	.size	_ZL6storeaU8__vectordPdi, .-_ZL6storeaU8__vectordPdi
	.type	_ZL6storeuU8__vectordPdi, @function
_ZL6storeuU8__vectordPdi:
.LFB856:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$12, %rsp
	vmovapd	%ymm0, -108(%rsp)
	movq	%rdi, -116(%rsp)
	movl	%esi, -120(%rsp)
	movl	-120(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-116(%rsp), %rax
	addq	%rdx, %rax
	movq	%rax, -20(%rsp)
	vmovapd	-108(%rsp), %ymm0
	vmovapd	%ymm0, -76(%rsp)
	movq	-20(%rsp), %rax
	vmovapd	-76(%rsp), %ymm0
	vmovupd	%ymm0, (%rax)
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE856:
	.size	_ZL6storeuU8__vectordPdi, .-_ZL6storeuU8__vectordPdi
	.type	_ZL5bcastPdi, @function
_ZL5bcastPdi:
.LFB858:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	movq	%rdi, -40(%rsp)
	movl	%esi, -44(%rsp)
	movl	-44(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rsp), %rax
	addq	%rdx, %rax
	movq	%rax, -24(%rsp)
	movq	-24(%rsp), %rax
	vbroadcastsd	(%rax), %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE858:
	.size	_ZL5bcastPdi, .-_ZL5bcastPdi
	.type	_ZL5bcastl, @function
_ZL5bcastl:
.LFB859:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	movq	%rdi, -40(%rsp)
	movq	-40(%rsp), %rax
	movq	%rax, -24(%rsp)
	vbroadcastsd	-24(%rsp), %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE859:
	.size	_ZL5bcastl, .-_ZL5bcastl
	.type	_ZL8mult_addU8__vectordS_S_, @function
_ZL8mult_addU8__vectordS_S_:
.LFB861:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	subq	$104, %rsp
	vmovapd	%ymm0, -56(%rsp)
	vmovapd	%ymm1, -88(%rsp)
	vmovapd	%ymm2, -120(%rsp)
	vmovapd	-56(%rsp), %ymm0
	vmovapd	%ymm0, 72(%rsp)
	vmovapd	-88(%rsp), %ymm0
	vmovapd	%ymm0, 40(%rsp)
	vmovapd	40(%rsp), %ymm0
	vmovapd	72(%rsp), %ymm1
	vmulpd	%ymm0, %ymm1, %ymm0
	vmovapd	%ymm0, 8(%rsp)
	vmovapd	-120(%rsp), %ymm0
	vmovapd	%ymm0, -24(%rsp)
	vmovapd	-24(%rsp), %ymm0
	vmovapd	8(%rsp), %ymm1
	vaddpd	%ymm0, %ymm1, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE861:
	.size	_ZL8mult_addU8__vectordS_S_, .-_ZL8mult_addU8__vectordS_S_
	.type	_ZL8mult_subU8__vectordS_S_, @function
_ZL8mult_subU8__vectordS_S_:
.LFB862:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	subq	$104, %rsp
	vmovapd	%ymm0, -56(%rsp)
	vmovapd	%ymm1, -88(%rsp)
	vmovapd	%ymm2, -120(%rsp)
	vmovapd	-56(%rsp), %ymm0
	vmovapd	%ymm0, 72(%rsp)
	vmovapd	-88(%rsp), %ymm0
	vmovapd	%ymm0, 40(%rsp)
	vmovapd	40(%rsp), %ymm0
	vmovapd	72(%rsp), %ymm1
	vmulpd	%ymm0, %ymm1, %ymm0
	vmovapd	-120(%rsp), %ymm1
	vmovapd	%ymm1, 8(%rsp)
	vmovapd	%ymm0, -24(%rsp)
	vmovapd	-24(%rsp), %ymm1
	vmovapd	8(%rsp), %ymm0
	vsubpd	%ymm1, %ymm0, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE862:
	.size	_ZL8mult_subU8__vectordS_S_, .-_ZL8mult_subU8__vectordS_S_
	.type	_ZL4multU8__vectordS_, @function
_ZL4multU8__vectordS_:
.LFB863:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	subq	$8, %rsp
	vmovapd	%ymm0, -88(%rsp)
	vmovapd	%ymm1, -120(%rsp)
	vmovapd	-88(%rsp), %ymm0
	vmovapd	%ymm0, -24(%rsp)
	vmovapd	-120(%rsp), %ymm0
	vmovapd	%ymm0, -56(%rsp)
	vmovapd	-56(%rsp), %ymm0
	vmovapd	-24(%rsp), %ymm1
	vmulpd	%ymm0, %ymm1, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE863:
	.size	_ZL4multU8__vectordS_, .-_ZL4multU8__vectordS_
	.type	_ZL3addU8__vectordS_, @function
_ZL3addU8__vectordS_:
.LFB864:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	subq	$8, %rsp
	vmovapd	%ymm0, -88(%rsp)
	vmovapd	%ymm1, -120(%rsp)
	vmovapd	-88(%rsp), %ymm0
	vmovapd	%ymm0, -24(%rsp)
	vmovapd	-120(%rsp), %ymm0
	vmovapd	%ymm0, -56(%rsp)
	vmovapd	-56(%rsp), %ymm0
	vmovapd	-24(%rsp), %ymm1
	vaddpd	%ymm0, %ymm1, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE864:
	.size	_ZL3addU8__vectordS_, .-_ZL3addU8__vectordS_
	.type	_ZL6cmpgtrU8__vectordS_, @function
_ZL6cmpgtrU8__vectordS_:
.LFB866:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	vmovapd	%ymm0, -48(%rsp)
	vmovapd	%ymm1, -80(%rsp)
	vmovapd	-48(%rsp), %ymm0
	vmovapd	-80(%rsp), %ymm1
	vcmppd	$14, %ymm1, %ymm0, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE866:
	.size	_ZL6cmpgtrU8__vectordS_, .-_ZL6cmpgtrU8__vectordS_
	.type	_ZL7select1U8__vectordS_, @function
_ZL7select1U8__vectordS_:
.LFB869:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	vmovapd	%ymm0, -48(%rsp)
	vmovapd	%ymm1, -80(%rsp)
	vmovapd	-48(%rsp), %ymm0
	vmovapd	-80(%rsp), %ymm1
	vblendpd	$1, %ymm1, %ymm0, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE869:
	.size	_ZL7select1U8__vectordS_, .-_ZL7select1U8__vectordS_
	.type	_ZL8select12U8__vectordS_, @function
_ZL8select12U8__vectordS_:
.LFB871:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	vmovapd	%ymm0, -48(%rsp)
	vmovapd	%ymm1, -80(%rsp)
	vmovapd	-48(%rsp), %ymm0
	vmovapd	-80(%rsp), %ymm1
	vblendpd	$3, %ymm1, %ymm0, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE871:
	.size	_ZL8select12U8__vectordS_, .-_ZL8select12U8__vectordS_
	.type	_ZL9select123U8__vectordS_, @function
_ZL9select123U8__vectordS_:
.LFB873:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	vmovapd	%ymm0, -48(%rsp)
	vmovapd	%ymm1, -80(%rsp)
	vmovapd	-48(%rsp), %ymm0
	vmovapd	-80(%rsp), %ymm1
	vblendpd	$7, %ymm1, %ymm0, %ymm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE873:
	.size	_ZL9select123U8__vectordS_, .-_ZL9select123U8__vectordS_
	.section	.rodata
.LC0:
	.string	"%lf %lf %lf %lf\n"
	.text
	.type	_ZL11printDArrayPd, @function
_ZL11printDArrayPd:
.LFB875:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rsi
	movq	-8(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rcx
	movq	-8(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	%rsi, -16(%rbp)
	vmovsd	-16(%rbp), %xmm3
	movq	%rcx, -16(%rbp)
	vmovsd	-16(%rbp), %xmm2
	movq	%rdx, -16(%rbp)
	vmovsd	-16(%rbp), %xmm1
	movq	%rax, -16(%rbp)
	vmovsd	-16(%rbp), %xmm0
	movl	$.LC0, %edi
	movl	$4, %eax
	call	printf
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE875:
	.size	_ZL11printDArrayPd, .-_ZL11printDArrayPd
	.section	.rodata
.LC1:
	.string	"%d %d %d %d\n"
	.text
	.type	_ZL13printIntArrayPi, @function
_ZL13printIntArrayPi:
.LFB876:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$12, %rax
	movl	(%rax), %esi
	movq	-8(%rbp), %rax
	addq	$8, %rax
	movl	(%rax), %ecx
	movq	-8(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-8(%rbp), %rax
	movl	(%rax), %eax
	movl	%esi, %r8d
	movl	%eax, %esi
	movl	$.LC1, %edi
	movl	$0, %eax
	call	printf
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE876:
	.size	_ZL13printIntArrayPi, .-_ZL13printIntArrayPi
	.type	_ZL10printm256IU8__vectorx, @function
_ZL10printm256IU8__vectorx:
.LFB877:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$-128, %rsp
	vmovdqa	%ymm0, (%rsp)
	leaq	32(%rsp), %rax
	movq	%rax, 120(%rsp)
	vmovdqa	(%rsp), %ymm0
	vmovdqa	%ymm0, 64(%rsp)
	movq	120(%rsp), %rax
	vmovdqa	64(%rsp), %ymm0
	vmovdqa	%ymm0, (%rax)
	movq	56(%rsp), %rsi
	movq	48(%rsp), %rcx
	movq	40(%rsp), %rdx
	movq	32(%rsp), %rax
	movq	%rsi, %r8
	movq	%rax, %rsi
	movl	$.LC1, %edi
	movl	$0, %eax
	call	printf
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE877:
	.size	_ZL10printm256IU8__vectorx, .-_ZL10printm256IU8__vectorx
	.type	_ZL10printm256DU8__vectord, @function
_ZL10printm256DU8__vectord:
.LFB878:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	subq	$160, %rsp
	vmovapd	%ymm0, 32(%rsp)
	leaq	64(%rsp), %rax
	movq	%rax, 152(%rsp)
	vmovapd	32(%rsp), %ymm0
	vmovapd	%ymm0, 96(%rsp)
	movq	152(%rsp), %rax
	vmovapd	96(%rsp), %ymm0
	vmovapd	%ymm0, (%rax)
	movq	88(%rsp), %rsi
	movq	80(%rsp), %rcx
	movq	72(%rsp), %rdx
	movq	64(%rsp), %rax
	movq	%rsi, 24(%rsp)
	vmovsd	24(%rsp), %xmm3
	movq	%rcx, 24(%rsp)
	vmovsd	24(%rsp), %xmm2
	movq	%rdx, 24(%rsp)
	vmovsd	24(%rsp), %xmm1
	movq	%rax, 24(%rsp)
	vmovsd	24(%rsp), %xmm0
	movl	$.LC0, %edi
	movl	$4, %eax
	call	printf
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE878:
	.size	_ZL10printm256DU8__vectord, .-_ZL10printm256DU8__vectord
	.type	_ZL16makeshiftmultaddPdS_S_, @function
_ZL16makeshiftmultaddPdS_S_:
.LFB879:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	movl	$0, -4(%rbp)
	jmp	.L43
.L44:
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rbp), %rax
	addq	%rdx, %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm1
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-32(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, (%rax)
	addl	$1, -4(%rbp)
.L43:
	cmpl	$2, -4(%rbp)
	jle	.L44
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE879:
	.size	_ZL16makeshiftmultaddPdS_S_, .-_ZL16makeshiftmultaddPdS_S_
	.type	_ZL16makeshiftmultsubPdS_S_, @function
_ZL16makeshiftmultsubPdS_S_:
.LFB880:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	movl	$0, -4(%rbp)
	jmp	.L46
.L47:
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rbp), %rax
	addq	%rdx, %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm1
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm2
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-32(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm0
	vmulsd	%xmm0, %xmm2, %xmm0
	vsubsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, (%rax)
	addl	$1, -4(%rbp)
.L46:
	cmpl	$2, -4(%rbp)
	jle	.L47
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE880:
	.size	_ZL16makeshiftmultsubPdS_S_, .-_ZL16makeshiftmultsubPdS_S_
	.type	_ZL13makeshiftmultPdS_S_, @function
_ZL13makeshiftmultPdS_S_:
.LFB881:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	movl	$0, -4(%rbp)
	jmp	.L49
.L50:
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rbp), %rax
	addq	%rdx, %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm1
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-32(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, (%rax)
	addl	$1, -4(%rbp)
.L49:
	cmpl	$2, -4(%rbp)
	jle	.L50
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE881:
	.size	_ZL13makeshiftmultPdS_S_, .-_ZL13makeshiftmultPdS_S_
	.type	_ZL12makeshiftaddPdS_S_, @function
_ZL12makeshiftaddPdS_S_:
.LFB882:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	movl	$0, -4(%rbp)
	jmp	.L52
.L53:
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rbp), %rax
	addq	%rdx, %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm1
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-32(%rbp), %rdx
	addq	%rcx, %rdx
	vmovsd	(%rdx), %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, (%rax)
	addl	$1, -4(%rbp)
.L52:
	cmpl	$2, -4(%rbp)
	jle	.L53
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE882:
	.size	_ZL12makeshiftaddPdS_S_, .-_ZL12makeshiftaddPdS_S_
	.type	_ZL16makeshiftcomparePdS_S_, @function
_ZL16makeshiftcomparePdS_S_:
.LFB883:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	movl	$0, -4(%rbp)
	jmp	.L55
.L59:
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rbp), %rax
	addq	%rax, %rdx
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rcx
	movq	-32(%rbp), %rax
	addq	%rcx, %rax
	vmovsd	(%rax), %xmm0
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rcx
	movq	-24(%rbp), %rax
	addq	%rcx, %rax
	vmovsd	(%rax), %xmm1
	vucomisd	%xmm1, %xmm0
	jbe	.L61
	movabsq	$4607182418800017408, %rax
	jmp	.L58
.L61:
	movl	$0, %eax
.L58:
	movq	%rax, (%rdx)
	addl	$1, -4(%rbp)
.L55:
	cmpl	$2, -4(%rbp)
	jle	.L59
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE883:
	.size	_ZL16makeshiftcomparePdS_S_, .-_ZL16makeshiftcomparePdS_S_
	.section	.rodata
.LC4:
	.string	"expected output"
	.text
	.type	_ZL6expectv, @function
_ZL6expectv:
.LFB884:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	$.LC4, %edi
	call	puts
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE884:
	.size	_ZL6expectv, .-_ZL6expectv
	.section	.rodata
.LC5:
	.string	"actual output"
	.text
	.type	_ZL6actualv, @function
_ZL6actualv:
.LFB885:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	$.LC5, %edi
	call	puts
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE885:
	.size	_ZL6actualv, .-_ZL6actualv
	.section	.rodata
.LC13:
	.string	"testing broadcast"
.LC14:
	.string	"expecting %lf\n"
.LC15:
	.string	"expecting %lu\n"
.LC16:
	.string	"comparing blends (swizzles)"
.LC17:
	.string	"array A: "
.LC18:
	.string	"array B: "
.LC19:
	.string	"select a 0 only"
.LC20:
	.string	"select a 0 and 1"
.LC21:
	.string	"select a 0 through 2"
	.text
	.globl	main
	.type	main, @function
main:
.LFB886:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	subq	$368, %rsp
	movabsq	$4607182418800017408, %rax
	movq	%rax, 16(%rsp)
	movabsq	$4617540697942969549, %rax
	movq	%rax, 24(%rsp)
	movabsq	$4613937818241073152, %rax
	movq	%rax, 32(%rsp)
	movabsq	$4619567317775286272, %rax
	movq	%rax, 40(%rsp)
	movabsq	$4608083138725491507, %rax
	movq	%rax, 48(%rsp)
	movabsq	$4611686018427387904, %rax
	movq	%rax, 56(%rsp)
	movabsq	$4621819117588971520, %rax
	movq	%rax, 64(%rsp)
	movabsq	$4622382067542392832, %rax
	movq	%rax, 72(%rsp)
	movq	$10, 80(%rsp)
	movq	$22, 88(%rsp)
	movq	$43, 96(%rsp)
	movq	$54, 104(%rsp)
	movq	$3, 112(%rsp)
	movq	$4, 120(%rsp)
	movq	$6, 128(%rsp)
	movq	$7, 136(%rsp)
	leaq	16(%rsp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZL5loadaPdi
	vmovapd	%ymm0, 336(%rsp)
	leaq	48(%rsp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZL5loaduPdi
	vmovapd	%ymm0, 304(%rsp)
	call	_ZL6expectv
	leaq	16(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	leaq	48(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	call	_ZL6actualv
	vmovapd	336(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	vmovapd	304(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	call	_ZL6expectv
	leaq	16(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	leaq	48(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	leaq	16(%rsp), %rax
	vmovapd	336(%rsp), %ymm0
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZL6storeaU8__vectordPdi
	leaq	48(%rsp), %rax
	vmovapd	304(%rsp), %ymm0
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZL6storeuU8__vectordPdi
	call	_ZL6actualv
	leaq	16(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	leaq	48(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	leaq	16(%rsp), %rax
	movl	$1, %esi
	movq	%rax, %rdi
	call	_ZL5bcastPdi
	vmovapd	%ymm0, 272(%rsp)
	movq	96(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL5bcastl
	vmovdqa	%ymm0, 240(%rsp)
	movl	$.LC13, %edi
	movl	$0, %eax
	call	printf
	movq	24(%rsp), %rax
	movq	%rax, 8(%rsp)
	vmovsd	8(%rsp), %xmm0
	movl	$.LC14, %edi
	movl	$1, %eax
	call	printf
	movq	96(%rsp), %rax
	movq	%rax, %rsi
	movl	$.LC15, %edi
	movl	$0, %eax
	call	printf
	vmovapd	272(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	vmovdqa	240(%rsp), %ymm0
	call	_ZL10printm256IU8__vectorx
	vmovapd	336(%rsp), %ymm2
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL8mult_addU8__vectordS_S_
	vmovapd	%ymm0, 208(%rsp)
	vmovapd	336(%rsp), %ymm2
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL8mult_subU8__vectordS_S_
	vmovapd	%ymm0, 176(%rsp)
	leaq	144(%rsp), %rdx
	leaq	48(%rsp), %rcx
	leaq	16(%rsp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZL16makeshiftmultaddPdS_S_
	call	_ZL6expectv
	leaq	144(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	leaq	144(%rsp), %rdx
	leaq	48(%rsp), %rcx
	leaq	16(%rsp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZL16makeshiftmultsubPdS_S_
	leaq	144(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	call	_ZL6actualv
	vmovapd	208(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	vmovapd	176(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	call	_ZL6expectv
	leaq	144(%rsp), %rdx
	leaq	48(%rsp), %rcx
	leaq	16(%rsp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZL13makeshiftmultPdS_S_
	call	_ZL6actualv
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL4multU8__vectordS_
	vmovapd	%ymm0, 208(%rsp)
	vmovapd	208(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	leaq	144(%rsp), %rdx
	leaq	48(%rsp), %rcx
	leaq	16(%rsp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZL12makeshiftaddPdS_S_
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL3addU8__vectordS_
	vmovapd	%ymm0, 208(%rsp)
	call	_ZL6expectv
	leaq	144(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	call	_ZL6actualv
	vmovapd	208(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	call	_ZL6expectv
	leaq	144(%rsp), %rdx
	leaq	48(%rsp), %rcx
	leaq	16(%rsp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZL16makeshiftcomparePdS_S_
	leaq	144(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL6cmpgtrU8__vectordS_
	vmovapd	%ymm0, 208(%rsp)
	call	_ZL6actualv
	vmovapd	208(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	call	_ZL6expectv
	movl	$.LC16, %edi
	call	puts
	movl	$.LC17, %edi
	movl	$0, %eax
	call	printf
	leaq	16(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	movl	$.LC18, %edi
	movl	$0, %eax
	call	printf
	leaq	48(%rsp), %rax
	movq	%rax, %rdi
	call	_ZL11printDArrayPd
	movl	$.LC19, %edi
	call	puts
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL7select1U8__vectordS_
	vmovapd	%ymm0, 208(%rsp)
	vmovapd	208(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	movl	$.LC20, %edi
	call	puts
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL8select12U8__vectordS_
	vmovapd	%ymm0, 208(%rsp)
	vmovapd	208(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	movl	$.LC21, %edi
	call	puts
	vmovapd	304(%rsp), %ymm1
	vmovapd	336(%rsp), %ymm0
	call	_ZL9select123U8__vectordS_
	vmovapd	%ymm0, 208(%rsp)
	vmovapd	208(%rsp), %ymm0
	call	_ZL10printm256DU8__vectord
	movl	$0, %eax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE886:
	.size	main, .-main
	.ident	"GCC: (GNU) 4.8.5 20150623 (Red Hat 4.8.5-16)"
	.section	.note.GNU-stack,"",@progbits

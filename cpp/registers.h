


// preserve required registers
.macro save_registers
stp x28, x27, [sp, -0x60]!              // store at sp - 0x60, sp modified, sp = sp - 0x60
stp x26, x25, [sp, 0x10]                // store at sp + 0x10, sp not modified
stp x24, x23, [sp, 0x20]
stp x22, x21, [sp, 0x30]
stp x20, x19, [sp, 0x40]
stp x29, x30, [sp, 0x50]
add x29, sp, 0x50                       // new frame pointer
.endm
// restore required registers
.macro restore_registers
ldp x29, x30, [sp, 0x50]
ldp x20, x19, [sp, 0x40]
ldp x22, x21, [sp, 0x30]
ldp x24, x23, [sp, 0x20]
ldp x26, x25, [sp, 0x10]
ldp x28, x27, [sp], 0x60                // load from sp, sp modified, sp = sp + 0x60
.endm



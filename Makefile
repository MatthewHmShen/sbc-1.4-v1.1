# 定义编译器
CC = gcc
CFLAGS = -I./sbc -Wall -g

# 定义目录
SBC_DIR = sbc
SRC_DIR = src
OBJ_DIR = obj

# 定义目标可执行文件（默认编译 sbc_enc 和 sbc_dec）
TARGETS = sbc_enc sbc_dec

# 如果需要编译 sbctester，取消下面这行的注释
# TARGETS += sbc_tester

# sbc 目录下的公共库文件
SBC_SRCS = $(SBC_DIR)/sbc.c $(SBC_DIR)/sbc_primitives.c
SBC_OBJS = $(patsubst $(SBC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SBC_SRCS))

# src 目录下的独立程序
ENC_SRC = $(SRC_DIR)/sbcenc.c
ENC_OBJ = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(ENC_SRC))

DEC_SRC = $(SRC_DIR)/sbcdec.c
DEC_OBJ = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(DEC_SRC))

# 如果需要 sbctester，取消下面两行的注释
# TESTER_SRC = $(SRC_DIR)/sbctester.c
# TESTER_OBJ = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(TESTER_SRC))

# 默认编译所有目标
all: $(TARGETS)

# 编译 sbc_enc（编码器）
sbc_enc: $(ENC_OBJ) $(SBC_OBJS)
	$(CC) $^ -o $@

# 编译 sbc_dec（解码器）
sbc_dec: $(DEC_OBJ) $(SBC_OBJS)
	$(CC) $^ -o $@

# 如果需要 sbctester，取消下面规则的注释
# sbc_tester: $(TESTER_OBJ) $(SBC_OBJS)
#	$(CC) $^ -o $@

# 编译 sbc 库文件
$(OBJ_DIR)/%.o: $(SBC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

# 编译 src 下的程序文件
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

# 清理
clean:
	rm -rf $(OBJ_DIR) $(TARGETS)

.PHONY: all clean


#ifndef __COMMON_DATATYPES_H__
#define __COMMON_DATATYPES_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef long long int64_t; /**< @brief Signed 64-bit integer. */
typedef int int32_t;        /**< @brief Signed 32-bit integer. */
typedef short int16_t;     /**< @brief Signed 16-bit integer. */
typedef char int8_t;       /**< @brief Signed  8-bit integer. */

typedef long long unsigned uint64_t; /**< @brief Unsigned 64-bit integer. */
typedef int unsigned uint32_t;       /**< @brief Unsigned 32-bit integer. */
typedef short unsigned uint16_t;     /**< @brief Unsigned 16-bit integer. */
typedef char unsigned uint8_t;       /**< @brief Unsigned  8-bit integer. */

typedef float float32_t;
typedef char u8_status;

#define nullptr 0x00

#ifdef __cplusplus
}
#endif

#endif
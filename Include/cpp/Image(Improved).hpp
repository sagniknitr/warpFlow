#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdint.h>
#include <vector>
 
#include "image.hpp"
#include "log.hpp"
#include "utils.hpp"
 
using namespace seva::utils;
using namespace seva::graph;
using namespace seva::system;
 
const char *Image::TAG = "OFI_Image";
 
static uint32_t GetChannelSize(Image::Format format) {
    switch (format) {
    case Image::Format::U16:
    case Image::Format::S16:
    case Image::Format::F16:
    case Image::Format::YUV422_1P:
    case Image::Format::YUV422_UV:
    case Image::Format::YUV420_NV21_UV:
    case Image::Format::YUV420_NV12_UV: return 2;
    case Image::Format::U32:
    case Image::Format::S32:
    case Image::Format::F32: return 4;
    default: return 1;
    }
}
 
static const char *ImageFormatToString(Image::Format format) {
    switch (format) {
    case Image::Format::U8: return "U8";
    case Image::Format::S8: return "S8";
    case Image::Format::F8: return "F8";
    case Image::Format::U16: return "U16";
    case Image::Format::S16: return "S16";
    case Image::Format::F16: return "F16";
    case Image::Format::U32: return "U32";
    case Image::Format::S32: return "S32";
    case Image::Format::F32: return "F32";
    case Image::Format::RGBX888: return "RGBX888";
    case Image::Format::BGRX888: return "BGRX888";
    case Image::Format::RGB888: return "RGB888";
    case Image::Format::BGR888: return "BGR888";
    case Image::Format::LUV888: return "LUV888";
    case Image::Format::YUV422_1P: return "YUV422_1P";
    case Image::Format::YUV422_2P: return "YUV422_2P";
    case Image::Format::YUV422_UV: return "YUV422_UV";
    case Image::Format::YUV420_NV21: return "YUV420_NV21";
    case Image::Format::YUV420_NV12: return "YUV420_NV12";
    case Image::Format::YUV420_NV21_UV: return "YUV420_NV21_UV";
    case Image::Format::YUV420_NV12_UV: return "YUV420_NV12_UV";
    default: Log::E("Image", "Unsupported Image format %d", format); return "UNDEFINED";
    }
}
 
static bool IsSupportedFileExtension(const std::string &extension) {
    if (extension.compare("PGM") == 0 || extension.compare("pgm") == 0 || extension.compare("RGB") == 0 ||
        extension.compare("rgb") == 0 || extension.compare("yuv")) {
        return true;
    }
 
    Log::E("Image", "%s Unsupported image file format", __FUNCTION__);
    return false;
}
 
/* Set image's mMetaData
 * channelSize: size of channel
 * channel: dims of channel for this plane
 * width: width for this plane
 * height: width for this plane */
void Image::SetBufferMetaData(uint32_t plane, size_t channelSize, size_t channel, size_t width, size_t height) {
    GetMetaData().GetStrides().resize(GetMetaData().GetNumPlanes(), std::vector<unsigned int>(SEVA_DIM_MAX, 0));
    GetMetaData().GetDimensions().resize(GetMetaData().GetNumPlanes(), std::vector<unsigned int>(SEVA_DIM_MAX, 0));
 
    GetMetaData().SetNumDimensions(SEVA_DIM_MAX);
    GetMetaData().SetDimension(plane, SEVA_DIM_C, channel);
    GetMetaData().SetDimension(plane, SEVA_DIM_X, width);
    GetMetaData().SetDimension(plane, SEVA_DIM_Y, height);
 
    GetMetaData().SetStride(plane, SEVA_DIM_C, channelSize);
    GetMetaData().SetStride(plane, SEVA_DIM_X, channel * channelSize);
    GetMetaData().SetStride(plane, SEVA_DIM_Y, channel * channelSize * width);
}
 
Image::Image() {
    mIsVirtual = true;
    mAccessOK = false;
    mType = Buffer::Type::IMAGE;
    mFormat = Image::Format::UNKNOWN;
}
 
Image::Image(Image::Format format) {
    mIsVirtual = true;
    mAccessOK = false;
    mType = Buffer::Type::IMAGE;
    mFormat = format;
}
 
bool Image::Initialize(uint32_t width, uint32_t height, Image::Format format) {
    uint32_t channelsize;
 
    mWidth = width;
    mHeight = height;
    mRealWidth = width;
    mRealHeight = height;
    mFormat = format;
 
    mValidRegion.SetStartx(0);
    mValidRegion.SetStarty(0);
    mValidRegion.SetEndx(width);
    mValidRegion.SetEndy(height);
 
    channelsize = GetChannelSize(format);
 
    switch (format) {
    case Format::U8:
    case Format::S8:
    case Format::F8:
        /* don't care channel parameter */
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        break;
    case Format::U16:
    case Format::S16:
    case Format::F16:
        /* don't care channel parameter */
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        break;
    case Format::U32:
    case Format::S32:
    case Format::F32:
        /* don't care channel parameter */
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        break;
    case Format::RGBX888:
    case Format::BGRX888:
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 4, mWidth, mHeight);
        break;
    case Format::RGB888:
    case Format::BGR888:
    case Format::LUV888:
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 3, mWidth, mHeight);
        break;
    case Format::YUV422_1P:
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        break;
    case Format::YUV422_2P:
        GetMetaData().SetNumPlanes(2);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        SetBufferMetaData(1, channelsize, 1, mWidth, mHeight);
        break;
    case Format::YUV422_UV:
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        break;
    case Format::YUV420_NV21:
    case Format::YUV420_NV12:
        GetMetaData().SetNumPlanes(2);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        SetBufferMetaData(1, channelsize, 2, mWidth / 2, mHeight / 2);
        break;
    case Format::YUV420_NV21_UV:
    case Format::YUV420_NV12_UV:
        GetMetaData().SetNumPlanes(1);
        SetBufferMetaData(0, channelsize, 1, mWidth, mHeight);
        break;
    default: Log::E(TAG, "Unsupported Image format %d", format); return false;
    }
    GetMetaData().SetNumDataInfos(GetMetaData().GetNumPlanes());
    GetMetaData().SetNumDimensions(3);
 
    AllocateMemory();
 
    return true;
}
 
Image::Image(uint32_t width, uint32_t height, Image::Format format) {
    mIsVirtual = false;
    mAccessOK = true;
    mType = Buffer::Type::IMAGE;
 
    Initialize(width, height, format);
}
 
Image::Image(int fd, void *addr, uint32_t width, uint32_t height, Image::Format format) {
    mFd = fd;
    mVa = addr;
    mIsVirtual = false;
    mAccessOK = true;
    mType = Buffer::Type::IMAGE;
 
    Initialize(width, height, format);
}
 
Image::Image(int fd,
             void *addr,
             uint32_t offset,
             uint32_t size,
             uint32_t width,
             uint32_t height,
             Image::Format format) {
    mSize = size;
    mOffset = offset;
    mFd = fd;
    mVa = addr;
    mIsVirtual = false;
    mAccessOK = true;
    mType = Buffer::Type::IMAGE;
 
    Initialize(width, height, format);
}
 
void Image::GetValidRegion(rect *region) {
    region->SetStartx(mValidRegion.GetStartx());
    region->SetStarty(mValidRegion.GetStarty());
    region->SetEndx(mValidRegion.GetEndx());
    region->SetEndy(mValidRegion.GetEndy());
}
 
bool Image::SetValidRegion(rect *region) {
    if (!mAccessOK) {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return false;
    }
 
    mValidRegion.SetStartx(region->GetStartx());
    mValidRegion.SetStarty(region->GetStarty());
    mValidRegion.SetEndx(region->GetEndx());
    mValidRegion.SetEndy(region->GetEndy());
    return true;
}
 
bool Image::ReadFromFile(const char *fileName) {
    if (fileName == nullptr) {
        Log::E(TAG, "%s(): fileName is nullptr", __FUNCTION__);
        return false;
    }
 
    Log::V(TAG, "%s(%s)", __FUNCTION__, fileName);
 
    uint32_t width, height;
    uint16_t maxpixvalue;
    unsigned int size = 0;
 
    if (!mAccessOK) {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return false;
    }
 
    std::size_t found = std::string(fileName).find_last_of(".");
    std::string extension = std::string(fileName).substr(found + 1);
 
    std::ifstream f;
    std::string line;
 
    f.open(fileName, std::ifstream::in);
    if (f.is_open() == false) {
        Log::E(TAG, "cannot open file %s", fileName);
        return false;
    }
 
    if (extension.compare("PGM") == 0 || extension.compare("pgm") == 0) {
        std::getline(f, line);  // PX
        std::getline(f, line);  // COMMENT
        std::getline(f, line);  // W, H
        std::sscanf(line.c_str(), "%u %u", &width, &height);
        std::getline(f, line);  // BPP
        std::sscanf(line.c_str(), "%hu ", &maxpixvalue);
 
        if (width > GetWidth() || height > GetHeight()) {
            Log::E(TAG, "Too big file. width or height is larger than Image object");
            f.close();
            return false;
        }
        if (maxpixvalue == 255) {
            if (GetFormat() != Image::Format::U8) {
                Log::E(TAG, "Unproper Image format, %s's format is U8", fileName);
                f.close();
                return false;
            }
        } else {
            if (GetFormat() != Image::Format::U16) {
                Log::E(TAG, "Unproper Image format, %s's format is U16", fileName);
                f.close();
                return false;
            }
        }
 
    } else if (extension.compare("rgb") == 0) {
        std::sscanf(fileName, "%*[^_]_%ux%u_%*s", &width, &height);
        if (width > GetWidth() || height > GetHeight()) {
            Log::E(TAG, "Too big file. width or height is larger than Image object");
            f.close();
            return false;
        }
    }
 
    for (unsigned int i = 0; i < GetMetaData().GetNumPlanes(); i++)
        size += GetAllocatedMemorySize(i);
 
    f.read((char *)mDataInfo->mVa, (int)size);
    f.close();
 
    return true;
}
 
bool Image::ReadFromFile(const char *fileName, Format format, int w, int h) {
    if (fileName == nullptr) {
        Log::E(TAG, "%s(): fileName is nullptr", __FUNCTION__);
        return false;
    }
 
    unsigned int size = 0;
    unsigned int offset = 0;
 
    if (!mAccessOK) {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return false;
    }
 
    std::ifstream f;
 
    f.open(fileName, std::ifstream::in);
    if (f.is_open() == false) {
        Log::E(TAG, "cannot open file %s", fileName);
        return false;
    }
 
    switch (format) {
    case Image::Format::U8:
    case Image::Format::S8:
    case Image::Format::F8:
        offset = 0;
        size = (unsigned int)(w * h);
        break;
    case Image::Format::U16:
    case Image::Format::S16:
    case Image::Format::F16:
        offset = 0;
        size = (unsigned int)(w * h * 2);
        break;
    case Image::Format::U32:
    case Image::Format::S32:
    case Image::Format::F32:
    case Image::Format::RGBX888:
    case Image::Format::BGRX888:
        offset = 0;
        size = (unsigned int)(w * h * 4);
        break;
    case Image::Format::RGB888:
    case Image::Format::BGR888:
    case Image::Format::LUV888:
        offset = 0;
        size = (unsigned int)(w * h * 3);
        break;
    case Image::Format::YUV422_1P:
    case Image::Format::YUV422_2P:
    case Image::Format::YUV420_NV21_UV:
    case Image::Format::YUV420_NV12_UV:
        offset = 0;
        size = (unsigned int)(w * h * 2);
        break;
    case Image::Format::YUV420_NV21:
    case Image::Format::YUV420_NV12:
        offset = 0;
        size = (unsigned int)(w * h + (w / 2) * h);
        break;
    default: Log::W(TAG, "%s not supported image format", __FUNCTION__); break;
    }
 
    // std::cout << "Format: " << ImageFormatToString(format) <<  "size: " << size << "offset: " << offset << std::endl;
    SetMetaData((unsigned int)w, (unsigned int)h);
    f.seekg(offset, f.beg);
    f.read((char *)mDataInfo->mVa, (int)size);
    f.close();
 
    return true;
}
 
bool Image::ReadFromFileWithMeta(const char *fileName) {
    if (fileName == nullptr) {
        Log::E(TAG, "%s(): fileName is nullptr", __FUNCTION__);
        return false;
    }
 
    uint32_t width, height;
    uint16_t maxpixvalue;
    Image::Format format;
    unsigned int size = 0;
 
    if (!mAccessOK) {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return false;
    }
 
    std::size_t found = std::string(fileName).find_last_of(".");
    std::string extension = std::string(fileName).substr(found + 1);
    if (IsSupportedFileExtension(extension) == false) {
        return false;
    }
 
    std::ifstream f;
    std::string line;
 
    f.open(fileName, std::ifstream::in);
    if (f.is_open() == false) {
        Log::E(TAG, "cannot open file %s", fileName);
        return false;
    }
 
    if (extension.compare("PGM") == 0 || extension.compare("pgm") == 0) {
        std::getline(f, line);  // PX
 
        std::getline(f, line);  // COMMENT
        std::getline(f, line);  // W, H
        std::sscanf(line.c_str(), "%u %u", &width, &height);
        std::getline(f, line);  // Maximum Pixel value
        std::sscanf(line.c_str(), "%hu ", &maxpixvalue);
 
        if (maxpixvalue == 255)
            format = Image::Format::U8;
        else
            format = Image::Format::U16;
 
        SetMetaData(width, height, format);
    } else if (extension.compare("rgb") == 0) {
        std::sscanf(fileName, "%*[^_]_%ux%u_%*s", &width, &height);
        format = Image::Format::RGB888;
        SetMetaData(width, height, format);
    }
 
    for (unsigned int i = 0; i < GetMetaData().GetNumPlanes(); i++)
        size += GetAllocatedMemorySize(i);
 
    f.read((char *)mDataInfo->mVa, (int)size);
    f.close();
 
    return true;
}
 
size_t Image::ReadFromBuffer(void *buffer, size_t size) {
    size_t image_buffer_size = 0;
 
    if (!mAccessOK) {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return 0;
    }
 
    image_buffer_size = GetAllocatedMemorySize();
 
    if (size > image_buffer_size) {
        Log::W(TAG, "%s size is bigger than image buffer size", __FUNCTION__);
        std::memcpy(GetBufferPtr(), buffer, image_buffer_size);
        return image_buffer_size;
    }
 
    std::memcpy(GetBufferPtr(), buffer, size);
    return size;
}
 
bool Image::WriteToFile(const char *fileName) {
    if (fileName == nullptr) {
        Log::E(TAG, "%s(): fileName is nullptr", __FUNCTION__);
        return false;
    }
 
    Log::V(TAG, "%s(%s)", __FUNCTION__, fileName);
 
    std::size_t found = std::string(fileName).find_last_of(".");
    std::string extension = std::string(fileName).substr(found + 1);
    unsigned int size = 0;
 
    if (!mAccessOK) {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return false;
    }
 
    if (Utils::CreateDirectory(fileName) == false) {
        Log::E(TAG, "%s(): Can't create directory(%s)", __FUNCTION__, fileName);
        return false;
    }
 
    std::FILE *outfile = std::fopen(fileName, "wb+");
 
    if (outfile == nullptr) {
        Log::E(TAG, "%s(): Can't create file(%s)", __FUNCTION__, fileName);
        return false;
    }
 
    if (extension.compare("PGM") == 0 || extension.compare("pgm") == 0) {
        std::fprintf(outfile, "P5\n# %s\n", fileName);
        std::fprintf(outfile, "%u %u\n", GetWidth(), GetHeight());
        if (GetFormat() == Format::U8) {
            std::fprintf(outfile, "255\n");
        } else if (GetFormat() == Format::U16) {
            std::fprintf(outfile, "65535\n");
        }
    }
 
    switch (GetFormat()) {
    case Image::Format::U8:
    case Image::Format::S8:
    case Image::Format::F8: size = (unsigned int)(mRealWidth * mRealHeight); break;
    case Image::Format::U16:
    case Image::Format::S16:
    case Image::Format::F16: size = (unsigned int)(mRealWidth * mRealHeight * 2); break;
    case Image::Format::U32:
    case Image::Format::S32:
    case Image::Format::F32:
    case Image::Format::RGBX888:
    case Image::Format::BGRX888: size = (unsigned int)(mRealWidth * mRealHeight * 4); break;
    case Image::Format::RGB888:
    case Image::Format::BGR888:
    case Image::Format::LUV888: size = (unsigned int)(mRealWidth * mRealHeight * 3); break;
    case Image::Format::YUV422_1P:
    case Image::Format::YUV422_2P:
    case Image::Format::YUV420_NV21_UV:
    case Image::Format::YUV420_NV12_UV: size = (unsigned int)(mRealWidth * mRealHeight * 2); break;
    case Image::Format::YUV420_NV21:
    case Image::Format::YUV420_NV12:
        size = (unsigned int)(mRealWidth * mRealHeight + (mRealWidth / 2) * mRealHeight);
        break;
    default: Log::W(TAG, "%s not supported image format", __FUNCTION__); break;
    }
 
    std::fwrite((char *)mDataInfo->mVa, sizeof(char), size, outfile);
    std::fflush(outfile);
    std::fclose(outfile);
 
    return true;
}
 
size_t Image::WriteToBuffer(void *ptr, size_t size) {
    size_t image_buffer_size = 0;
 
    if (!mAccessOK) {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return 0;
    }
 
    image_buffer_size = GetAllocatedMemorySize();
 
    if (size > image_buffer_size) {
        Log::W(TAG, "%s size is bigger than image buffer size", __FUNCTION__);
        std::memcpy(ptr, GetBufferPtr(), image_buffer_size);
        return image_buffer_size;
    }
 
    std::memcpy(ptr, GetBufferPtr(), size);
    return size;
}
 
std::shared_ptr<Image> Image::CloneMeta(Image &source) {
    auto item = std::make_shared<Image>(source.GetWidth(), source.GetHeight(), source.GetFormat());
    item->mValidRegion = source.mValidRegion;
    return item;
}
 
void Image::DumpImageInfo(void) {
    Log::V(TAG, "%s()", __FUNCTION__);
 
    if (mIsVirtual == true) {
        Log::V(TAG, "%s(): This image is a virtual image.", __FUNCTION__);
        return;
    }
 
    uint32_t planes, dims;
 
    Log::V(TAG, "====================Dump Image Information====================");
    Log::V(TAG, "Image format: %s", ImageFormatToString(mFormat));
 
    for (planes = 0; planes < GetMetaData().GetNumPlanes(); planes++) {
        for (dims = 0; dims < GetMetaData().GetNumDimensions(); dims++) {
            Log::V(TAG, "mDimensions[%d][%d]: %d", planes, dims, GetMetaData().GetDimensions()[planes][dims]);
            Log::V(TAG, "mStrides[%d][%d]: %d", planes, dims, GetMetaData().GetStrides()[planes][dims]);
        }
    }
    Log::V(TAG,
           "mValidRegion : (%d, %d) - (%d, %d)",
           mValidRegion.GetStartx(),
           mValidRegion.GetStarty(),
           mValidRegion.GetEndx(),
           mValidRegion.GetEndy());
    Log::V(TAG, "Virtual : [%c]", mIsVirtual ? 'Y' : 'N');
    Log::V(TAG, "Accessok : [%c]", mAccessOK ? 'Y' : 'N');
    Log::V(TAG, "====================Dump Image Information====================");
}
 
Image::~Image() {}
 
bool Image::SetMetaData(uint32_t width, uint32_t height, Format format) {
    if (mAccessOK) {
        Initialize(width, height, format);
        mChanged = true;
        return true;
    } else {
        Log::W(TAG, "%s try to access virtual image", __FUNCTION__);
        return false;
    }
}
 
bool Image::SetMetaData(uint32_t width, uint32_t height) {
    mRealWidth = width;
    mRealHeight = height;
    return true;
}
 
void *Image::GetImageBuffer() { return GetBufferPtr(); }
 
void *Image::GetImagePlaneBuffer(int plane) { return GetPlaneBufferPtr((unsigned int)plane); }
 
std::shared_ptr<Image> Image::MakeImage(uint32_t width, uint32_t height, Format format) {
    return std::make_shared<Image>(width, height, format);
}
 
std::shared_ptr<Image> Image::MakeImage(int fd, void *addr, uint32_t width, uint32_t height, Format format) {
    return std::make_shared<Image>(fd, addr, width, height, format);
}
 
std::shared_ptr<Image>
Image::MakeImage(int fd, void *addr, uint32_t offset, uint32_t size, uint32_t width, uint32_t height, Format format) {
    return std::make_shared<Image>(fd, addr, offset, size, width, height, format);
}
 
std::shared_ptr<Image> Image::MakeImage(Format format) { return std::make_shared<Image>(format); }
 
std::shared_ptr<Image> Image::MakeImage() { return std::make_shared<Image>(); }
 



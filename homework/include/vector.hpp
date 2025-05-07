#pragma once

#include "core.hpp"


namespace CMF
{

class Vec3
{
public:
    real_t x = 0.0, y = 0.0, z = 0.0;

    constexpr Vec3() noexcept = default;
    constexpr Vec3(const Vec3&) noexcept = default;
    constexpr Vec3(Vec3&&) noexcept = default;
    constexpr Vec3& operator=(const Vec3&) noexcept = default;
    constexpr Vec3& operator=(Vec3&&) noexcept = default;

    constexpr Vec3(real_t x, real_t y, real_t z) noexcept
        : x(x), y(y), z(z)
    {
    }

    static inline constexpr Vec3 null() noexcept { return Vec3(0.0, 0.0, 0.0); }

    inline constexpr Vec3 operator+(const Vec3& other) const noexcept
    {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    inline constexpr Vec3 operator-(const Vec3& other) const noexcept
    {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    inline constexpr Vec3 operator*(real_t scalar) const noexcept
    {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    inline constexpr Vec3 operator/(real_t scalar) const noexcept
    {
        return Vec3(x / scalar, y / scalar, z / scalar);
    }

    inline constexpr Vec3& operator+=(const Vec3& other) noexcept
    {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    constexpr Vec3& operator-=(const Vec3& other) noexcept
    {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    constexpr Vec3& operator*=(real_t scalar) noexcept
    {
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    constexpr Vec3& operator/=(real_t scalar) noexcept
    {
        x /= scalar; y /= scalar; z /= scalar;
        return *this;
    }

    [[nodiscard]] inline constexpr real_t dot(const Vec3& other) const noexcept
    {
        return x * other.x + y * other.y + z * other.z;
    }

    [[nodiscard]] constexpr Vec3 cross(const Vec3& other) const noexcept
    {
        return Vec3(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    [[nodiscard]] inline constexpr real_t magnitude() const noexcept
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    [[nodiscard]] constexpr Vec3 normalized() const noexcept
    {
        return *this / magnitude();
    }

};  // class Vec3
}  // namespace CMF

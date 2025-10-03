# constants.py
import numpy as np

# Физические константы в CGS
M_SUN = 1.989e33           # г
PC = 3.086e18              # см
YR = 3.156e7               # с
KM = 1e5                   # см
ERG = 1.0                  # эрг
EV = 1.602e-12             # эрг

# Астрофизические константы
GAMMA = 5/3                # Показатель адиабаты
MU = 0.6                   # Средний молекулярный вес
MP = 1.673e-24             # г (масса протона)
KB = 1.381e-16             # эрг/K
C = 2.998e10               # см/с
SIGMA_SB = 5.67e-5         # эрг/cm²/K⁴/s

class Normalization:
    def __init__(self, length_norm=PC, density_norm=1e-24, velocity_norm=1e8):
        self.length = length_norm
        self.density = density_norm
        self.velocity = velocity_norm
        self.time = length_norm / velocity_norm
        self.pressure = density_norm * velocity_norm**2
        self.energy = density_norm * velocity_norm**2 * length_norm**3
        
        # Проверка на разумные значения
        self._validate_norms()
    
    def _validate_norms(self):
        """Проверка физической осмысленности нормировок"""
        if self.length <= 0 or self.density <= 0 or self.velocity <= 0:
            raise ValueError("Нормировки должны быть положительными")
        
        # Проверка на переполнение
        if self.energy > 1e60 or self.energy < 1e20:
            print(f"Предупреждение: энергетическая нормировка {self.energy:.2e} эрг может быть нефизичной")

def dimensional_to_dimensionless(val, norm):
    """Перевод размерных величин в безразмерные с проверкой"""
    if hasattr(val, '__len__'):
        val = np.asarray(val)
        if np.any(val < 0):
            print("Предупреждение: отрицательные значения при нормировке")
        return val / norm
    else:
        if val < 0:
            print("Предупреждение: отрицательное значение при нормировке")
        return val / norm

def dimensionless_to_dimensional(val, norm):
    """Перевод безразмерных величин в размерные с проверкой"""
    if hasattr(val, '__len__'):
        val = np.asarray(val)
        result = val * norm
        if np.any(np.isnan(result)) or np.any(np.isinf(result)):
            print("Ошибка: NaN или Inf при преобразовании в размерные величины")
        return result
    else:
        result = val * norm
        if np.isnan(result) or np.isinf(result):
            print("Ошибка: NaN или Inf при преобразовании в размерные величины")
        return result

if __name__ == "__main__":
    print("Тестирование constants.py:")
    print(f"M_SUN = {M_SUN:.2e} г")
    print(f"PC = {PC:.2e} см") 
    print(f"YR = {YR:.2e} с")
    print(f"GAMMA = {GAMMA}")
    
    # Тест нормировок
    norm = Normalization()
    test_length = 2 * PC
    dimless = dimensional_to_dimensionless(test_length, norm.length)
    dim_back = dimensionless_to_dimensional(dimless, norm.length)
    print(f"Тест нормировки: {test_length:.2e} см -> {dimless:.2f} -> {dim_back:.2e} см")
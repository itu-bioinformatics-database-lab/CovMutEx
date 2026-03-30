import axios from "axios";

// Backend URL (Django)
const API_BASE_URL = "http://127.0.0.1:8000/api";

export const api = {
  predict: (data) => {
    // Eğer dosya yükleniyorsa (FormData)
    if (data.isNewUpload) {
      const formData = new FormData();
      formData.append("nodeId", data.nodeId);
      formData.append("elapsedDay", data.elapsedDay);
      formData.append("uploadFolderName", data.uploadName); // Backend bu ismi bekliyor
      if (data.selectedProteinRegion) {
        formData.append("selectedProteinRegion", data.selectedProteinRegion);
      }
      
      // Dosyalar
      if (data.modelFile) formData.append("modelFile", data.modelFile);
      if (data.extractorFile) formData.append("extractorFile", data.extractorFile);
      
      // Helper Files
      if (data.helperFiles && data.helperFiles.length > 0) {
        data.helperFiles.forEach((file, index) => {
          formData.append(`helperFile_${index}`, file);
          formData.append(`helperFileName_${index}`, file.name);
        });
      }

      // Custom Params (JSON String olarak gönderelim)
      if (data.customParameters) {
        formData.append("customParameters", JSON.stringify(data.customParameters));
      }

      return axios.post(`${API_BASE_URL}/predict/`, formData, {
        headers: { "Content-Type": "multipart/form-data" },
      });
    } 
    // Standart JSON isteği
    else {
      return axios.post(`${API_BASE_URL}/predict/`, data);
    }
  },

  getModels: () => axios.get(`${API_BASE_URL}/models/`),
  getModelParams: (modelName) => axios.get(`${API_BASE_URL}/model-parameters/`, { params: { model_name: modelName } }),
};
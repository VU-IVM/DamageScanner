// set-body-class.js
document.addEventListener("DOMContentLoaded", function () {
    const path = window.location.pathname.toLowerCase();
  
    // Match pages under "how-to/" or "how-it-works/"
    if (path.includes("/how-to/") || path.includes("/how-it-works/")) {
      document.body.classList.add("show-sidebar");
    }
  });
// set-body-class.js
document.addEventListener("DOMContentLoaded", function () {
  const path = window.location.pathname.toLowerCase();

  // Match pages under "how-to/", "how-it-works/", or "examples/"
  if (
    path.includes("/how-to/") ||
    path.includes("/how-it-works/") ||
    path.includes("/examples/")
  ) {
    document.body.classList.add("show-sidebar");
  }
});